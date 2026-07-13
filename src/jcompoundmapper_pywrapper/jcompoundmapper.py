"""Python wrapper for jCompoundMapper fingerprints"""

import logging
import multiprocessing
import os
import subprocess
import warnings
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from dataclasses import dataclass
from enum import Enum, auto

import more_itertools
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.rdBase import BlockLogs

from .utils import install_java, mktempfile, needsHs, read_libsvmsparse

logger = logging.getLogger(__name__)


class AtomType(Enum):
    CDK_ATOM_TYPES = auto()
    ELEMENT_NEIGHBOR = auto()
    ELEMENT_NEIGHBOR_RING = auto()
    CUSTOM = auto()
    ELEMENT_SYMBOL = auto()
    DAYLIGHT_INVARIANT = auto()
    DAYLIGHT_INVARIANT_RING = auto()


@dataclass
class FpParamType:
    """Parameters of a jCompoundMapper fingerprint.

    :param depth: only for 2D fingerprints: DFS, ASP, AP2D, AT2D, CATS2D, PHAP2POINT2D, PHAP3POINT2D,
    ECFP, ECFPVariant, LSTAR, SHED, RAD2D, MACCS
    :param dist_cutoff: only for 3D fingerprints: AP3D, AT3D, CATS3D, PHAP2POINT3D, PHAP3POINT3D, RAD3D
    :param stretch_factor: only for 3D fingerprints
    :param atom_type:
    :param arom_flag: default: False
    """

    depth: int | None
    dist_cutoff: float | None
    stretch_factor: float | None
    atom_type: AtomType
    arom_flag: bool


DEFAULT_FP_PARAMETERS = {'DFS': FpParamType(depth=8, dist_cutoff=None, stretch_factor=None,
                                            atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'ASP': FpParamType(depth=8, dist_cutoff=None, stretch_factor=None,
                                            atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'AP2D': FpParamType(depth=8, dist_cutoff=None, stretch_factor=None,
                                             atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'AT2D': FpParamType(depth=5, dist_cutoff=None, stretch_factor=None,
                                             atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'AP3D': FpParamType(depth=None, dist_cutoff=10, stretch_factor=1,
                                             atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'AT3D': FpParamType(depth=None, dist_cutoff=6, stretch_factor=1,
                                             atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'CATS2D': FpParamType(depth=9, dist_cutoff=None, stretch_factor=None,
                                               atom_type=AtomType.CUSTOM, arom_flag=False),
                         'CATS3D': FpParamType(depth=None, dist_cutoff=9, stretch_factor=1,
                                               atom_type=AtomType.CUSTOM, arom_flag=False),
                         'PHAP2POINT2D': FpParamType(depth=8, dist_cutoff=None, stretch_factor=None,
                                                     atom_type=AtomType.CUSTOM, arom_flag=False),
                         'PHAP3POINT2D': FpParamType(depth=5, dist_cutoff=None, stretch_factor=None,
                                                     atom_type=AtomType.CUSTOM, arom_flag=False),
                         'PHAP2POINT3D': FpParamType(depth=None, dist_cutoff=10, stretch_factor=1,
                                                     atom_type=AtomType.CUSTOM, arom_flag=False),
                         'PHAP3POINT3D': FpParamType(depth=None, dist_cutoff=6, stretch_factor=1,
                                                     atom_type=AtomType.CUSTOM, arom_flag=False),
                         'ECFP': FpParamType(depth=8, dist_cutoff=None, stretch_factor=None,
                                             atom_type=AtomType.DAYLIGHT_INVARIANT_RING,
                                             arom_flag=False),
                         'ECFPVariant': FpParamType(depth=4, dist_cutoff=None, stretch_factor=None,
                                                    atom_type=AtomType.DAYLIGHT_INVARIANT_RING,
                                                    arom_flag=False),
                         'LSTAR': FpParamType(depth=6, dist_cutoff=None, stretch_factor=None,
                                              atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'SHED': FpParamType(depth=8, dist_cutoff=None, stretch_factor=None,
                                             atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'RAD2D': FpParamType(depth=3, dist_cutoff=None, stretch_factor=None,
                                              atom_type=AtomType.ELEMENT_SYMBOL, arom_flag=False),
                         'RAD3D': FpParamType(depth=None, dist_cutoff=4, stretch_factor=1,
                                              atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False),
                         'MACCS': FpParamType(depth=8, dist_cutoff=None, stretch_factor=None,
                                              atom_type=AtomType.ELEMENT_NEIGHBOR, arom_flag=False)
                         }



Fingerprint = Enum('Fingerprint', list(DEFAULT_FP_PARAMETERS))

Fingerprints_2D = [fp_type for fp_type in Fingerprint if '3d' not in str(fp_type).lower()]
Fingerprints_3D = [fp_type for fp_type in Fingerprint if '3d' in str(fp_type).lower()]


def _make_chunks(mols: list, njobs: int, chunksize: int | None) -> list[list]:
    """Split molecules into chunks to be dispatched to worker processes.

    If chunksize is None, mols are split into exactly min(njobs, len(mols))
    balanced chunks (sizes differ by at most one) so every requested worker
    is guaranteed to receive work. If chunksize is given, it is authoritative
    and mols are batched into fixed-size chunks (may leave some workers idle).

    :param mols: molecules (or any list) to split into chunks
    :param njobs: number of concurrent workers that will consume the chunks
    :param chunksize: fixed size of each chunk, or None to auto-balance across njobs
    :return: a list of chunks
    """
    if not mols:
        return []
    if chunksize is None:
        n_chunks = min(njobs, len(mols))
        base, remainder = divmod(len(mols), n_chunks)
        chunks = []
        start = 0
        for i in range(n_chunks):
            size = base + (1 if i < remainder else 0)
            chunks.append(mols[start:start + size])
            start += size
        return chunks
    return [list(chunk) for chunk in more_itertools.batched(mols, chunksize)]


class JCompoundMapper:
    """Wrapper to obtain molecular fingerprints from jCompoundMapper.

    jCompoundMapper: An open source Java library and command-line tool for chemical fingerprints.
    Hinselmann, G., Rosenbaum, L., Jahn, A., Fechner N., and Zeel A., J Cheminform 3, 3 (2011).
    https://doi.org/10.1186/1758-2946-3-3
    """

    # Extra safety net for JRE installation; the real guarantee is that calculate() installs
    # the JRE once, before any worker starts.
    lock = multiprocessing.RLock()
    # Path to the JAR file
    _jarfile = os.path.abspath(os.path.join(__file__, os.pardir, 'jCompoundMapper', 'jCMapperCLI.jar'))

    def __init__(self, fingerprint: str | Fingerprint = 'DFS', params: FpParamType | None = None):
        """Instantiate a wrapper to calculate jCompoundMapper molecular fingerprints."""
        if not isinstance(fingerprint, Fingerprint) and fingerprint not in DEFAULT_FP_PARAMETERS.keys():
            raise ValueError('fingerprint can only be one of {' + ', '.join(DEFAULT_FP_PARAMETERS.keys()) + '}')
        elif params is not None and not isinstance(params, FpParamType):
            raise ValueError('fingerprint params must be an instance of FpParamType')
        # Ensure the jar file exists
        if not os.path.isfile(self._jarfile):
            raise OSError('The required JAR file is not present. Reinstall jcompoundmapper-pywrapper.')
        # Define internal parameters
        self.fp_name = fingerprint if isinstance(fingerprint, str) else fingerprint.name
        self.fp_params = params

    def calculate(self, mols: Iterable[Chem.Mol], nbits: int = 1024, show_banner: bool = True, njobs: int = 1,
                  chunksize: int | None = None) -> pd.DataFrame:
        """Calculate molecular fingerprints.

        :param mols: RDKit molecules for which fingerprints should be calculated
        :param nbits: size of fingerprints
        :param show_banner: If True, show notice on fingerprint usage
        :param njobs: number of concurrent processes; must not exceed the number of available CPU cores
        :param chunksize: number of molecules processed per worker process; if None (default) molecules
            are auto-balanced across njobs workers so every worker gets work; ignored if njobs is 1
        :return: a pandas DataFrame containing the fingerprint values
        """
        if njobs < 1:
            raise ValueError('njobs must be a strictly positive integer.')
        cpu_count = os.cpu_count() or 1
        if njobs > cpu_count:
            raise ValueError(f'njobs ({njobs}) exceeds the number of available CPU cores ({cpu_count}).')
        if self.fp_name == 'MACCS':
            nbits = 166
        if show_banner:
            self._show_banner()
        mols = list(mols)
        # Parallelize should need be
        if njobs > 1:
            # Ensure the JRE is installed once in the parent process before workers start,
            # so concurrently-spawned workers all find it already in place.
            with self.lock:
                install_java()
            chunks = _make_chunks(mols, njobs, chunksize)
            if not chunks:
                return pd.DataFrame()
            with ProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._multiproc_calculate, chunk, nbits) for chunk in chunks]
                results = [future.result() for future in futures]
            return pd.concat(results).reset_index(drop=True)
        # Single process
        return self._calculate(mols, nbits=nbits)

    def _show_banner(self):
        """Log info message for citing."""
        logger.info("""Fingerprinting algorithms for your use.
jCompoundMapper provides popular fingerprinting algorithms for chemical graphs such as
depth-first search fingerprints, shortest-path fingerprints, extended connectivity fingerprints,
autocorrelation fingerprints (e.g. CATS2D), radial fingerprints (e.g. Molprint2D), geometrical
Molprint, atom pairs, and pharmacophore fingerprints.

###################################

Should you publish results based on the jCompoundMapper fingerprints, please cite:

jCompoundMapper: An open source Java library and command-line tool for chemical fingerprints.
Hinselmann, G., Rosenbaum, L., Jahn, A., Fechner N., and Zeel A.
Journal of Cheminformatics 2003  3 (3).
DOI: 10.1186/1758-2946-3-3

###################################
""")

    def _prepare_command(self, mols: list[Chem.Mol], nbits: int = 1024) -> list[str]:
        """Create the JCompoundMapper command to be run to obtain molecular fingerprints.

        :param mols: molecules to obtained molecular fingerprint of
        :param nbits: the size of the output fingerprints
        :return: the command, as an argument list, to run.
        """
        # 1) Ensure JRE is accessible
        with self.lock:
            self._java_path = install_java()
        # 2) Create temp SD v2k file
        self._tmp_sd = mktempfile('molecules_v2k.sd')
        self._out = mktempfile('jcm_output.tsv')
        self._skipped = []
        try:
            with BlockLogs():
                writer = Chem.SDWriter(self._tmp_sd)
                # Ensure V2000 as jCompoundMapper cannot properly process v3000
                writer.SetForceV3000(False)
                for i, mol in enumerate(mols):
                    if mol is not None and isinstance(mol, Chem.Mol):
                        if mol.GetNumAtoms() > 999:
                            raise ValueError('Cannot calculate fingerprint for molecules with more than 999 atoms.')
                        # Does molecule lack hydrogen atoms?
                        if needsHs(mol):
                            warnings.warn(
                                'Molecule lacks hydrogen atoms: this might affect the value of calculated fingerprint',
                                category=UserWarning, stacklevel=2,
                            )
                        # 3D fingerprint
                        if '3D' in self.fp_name:
                            confs = list(mol.GetConformers())
                            if not (len(confs) > 0 and confs[-1].Is3D()):
                                raise ValueError('Cannot calculate the 3D fingerprint of a conformer-less molecule')
                        writer.write(mol)
                    else:
                        self._skipped.append(i)
                writer.close()
        except ValueError as e:
            # Free resources and raise error
            writer.close()
            os.remove(self._tmp_sd)
            raise e from None
        # 3) Create command
        command_params: list[str] = []
        if self.fp_params is not None:
            if self.fp_params.atom_type is not None:
                command_params += ['-a', self.fp_params.atom_type.name]
            if self.fp_params.depth is not None:
                command_params += ['-d', str(self.fp_params.depth)]
            if self.fp_params.dist_cutoff is not None:
                command_params += ['-d', str(self.fp_params.dist_cutoff)]
            if self.fp_params.arom_flag is not None:
                command_params += ['-k', str(self.fp_params.arom_flag).lower()]
            if self.fp_params.stretch_factor is not None:
                command_params += ['-s', str(self.fp_params.stretch_factor)]
        command = [
            self._java_path,
            '-Djava.awt.headless=true',
            # Pin each JVM to a single core: parallelism comes from spawning njobs OS
            # processes, not from letting each JVM auto-size GC/JIT threads to the host's
            # full core count (which would oversubscribe cores beyond what njobs requested).
            '-XX:ActiveProcessorCount=1',
            '-jar', self._jarfile,
            '-f', self._tmp_sd,
            '-ff', 'LIBSVM_SPARSE',
            '-hs', str(nbits),
            *command_params,
            '-o', self._out,
        ]
        return command

    def _cleanup(self, sd: bool = True, output: bool = True) -> None:
        """Cleanup resources used for calculation."""
        # Remove temporary files
        if sd:
            os.remove(self._tmp_sd)
        if output:
            os.remove(self._out)

    def _run_command(self, command: list[str], nbits: int = 1024) -> pd.DataFrame:
        """Run the jCompoundMapper command.

        :param command: the command, as an argument list, to be run
        :param nbits: size of output fingerprints
        """
        process = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if process.returncode == 0:
            values = read_libsvmsparse(self._out, nbits)
        else:
            stderr_text = process.stderr.decode(errors='replace') if process.stderr else ''
            self._cleanup()
            raise ValueError(
                f'A conflicting set of parameters was used. Please fix or use the defaults.\n{stderr_text}'
            )
        values = pd.DataFrame(values, columns=[f'{self.fp_name}_{i}' for i in range(1, nbits + 1)])
        return values

    def _calculate(self, mols: list[Chem.Mol], nbits: int = 1024) -> pd.DataFrame:
        """Calculate JCompoundMapper fingerprints on one process.

        :param mols: RDkit molecules for which JCompoundMapper fingerprints should be calculated.
        :param nbits: size of the fingerprint
        :return: a pandas DataFrame containing fingerprint values
        """
        # Prepare inputs
        command = self._prepare_command(mols, nbits)
        # Run command and obtain results
        results = self._run_command(command, nbits)
        # Cleanup
        self._cleanup()
        # Insert lines of skipped molecules
        if len(self._skipped):
            results = (pd.DataFrame(np.insert(results.values, self._skipped,
                                              values=[np.nan] * len(results.columns),
                                              axis=0),
                                    columns=results.columns)
                       )
        results = (results.apply(pd.to_numeric, errors='coerce')
                          .fillna(0)
                          .convert_dtypes()
                   )
        return results

    def _multiproc_calculate(self, mols: list[Chem.Mol], nbits: int = 1024) -> pd.DataFrame:
        """Calculate jCompoundMapper fingerprints in a separate worker process.

        :param mols: RDKit molecules for which jCompoundMapper fingerprints should be calculated
        :param nbits: size of the fingerprint
        :return: a pandas DataFrame containing the jCompoundMapper fingerprint values
        """
        # Copy self instance to make thread safe
        jcm = deepcopy(self)
        # Run copy
        result = jcm.calculate(mols, nbits=nbits, show_banner=False, njobs=1)
        return result
