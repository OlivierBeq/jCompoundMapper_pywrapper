# -*- coding: utf-8

"""Python wrapper for jCompoundMapper fingerprints"""

from __future__ import annotations

import os
import multiprocessing
import warnings
import subprocess
from copy import deepcopy
from enum import Enum, auto
from typing import Iterable, List, NamedTuple, Optional

import more_itertools
import numpy as np
import pandas as pd
from bounded_pool_executor import BoundedProcessPoolExecutor
from rdkit import Chem
from rdkit.rdBase import BlockLogs

from .utils import install_java, mktempfile, needsHs, read_libsvmsparse


class AtomType(Enum):
    CDK_ATOM_TYPES = auto()
    ELEMENT_NEIGHBOR = auto()
    ELEMENT_NEIGHBOR_RING = auto()
    CUSTOM = auto()
    ELEMENT_SYMBOL = auto()
    DAYLIGHT_INVARIANT = auto()
    DAYLIGHT_INVARIANT_RING = auto()


class FpParamType:

    def __init__(self, depth: Optional[int], dist_cutoff: Optional[float],
                 stretch_factor: Optional[float], atom_type: AtomType,
                 arom_flag: bool):
        """

        :param depth: only for 2D fingerprints: DFS, ASP, AP2D, AT2D, CATS2D, PHAP2POINT2D, PHAP3POINT2D,
        ECFP, ECFPVariant, LSTAR, SHED, RAD2D, MACCS
        :param dist_cutoff: only for 3D fingerprints: AP3D, AT3D, CATS3D, PHAP2POINT3D, PHAP3POINT3D, RAD3D
        :param stretch_factor: only for 3D fingerprints
        :param atom_type:
        :param arom_flag: default: False
        """
        self.depth = depth
        self.dist_cutoff = dist_cutoff
        self.stretch_factor = stretch_factor
        self.atom_type = atom_type
        self.arom_flag = arom_flag


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



Fingerprint = Enum('Fingerprint', [fp_type for fp_type in DEFAULT_FP_PARAMETERS.keys()])

Fingerprints_2D = [fp_type for fp_type in Fingerprint if '3d' not in str(fp_type).lower()]
Fingerprints_3D = [fp_type for fp_type in Fingerprint if '3d' in str(fp_type).lower()]

class JCompoundMapper:
    """Wrapper to obtain molecular fingerprints from jCompoundMapper.

    jCompoundMapper: An open source Java library and command-line tool for chemical fingerprints.
    Hinselmann, G., Rosenbaum, L., Jahn, A., Fechner N., and Zeel A., J Cheminform 3, 3 (2011).
    https://doi.org/10.1186/1758-2946-3-3
    """

    lock = multiprocessing.RLock()  # Ensure installation of JRE is thread safe
    _jarfile = os.path.abspath(os.path.join(__file__, os.pardir, 'jCompoundMapper', 'jCMapperCLI.jar'))  # Path to the JAR file

    def __init__(self, fingerprint: str | Fingerprint = 'DFS', params: Optional[FpParamType] = None):
        """Instantiate a wrapper to calculate jCompoundMapper molecular fingerprints."""
        if not isinstance(fingerprint, Fingerprint) and fingerprint not in DEFAULT_FP_PARAMETERS.keys():
            raise ValueError('fingerprint can only be one of {' + ', '.join(DEFAULT_FP_PARAMETERS.keys()) + '}')
        elif params is not None and not isinstance(params, FpParamType):
            raise ValueError('fingerprint params must be an instance of FpParamType')
        # Ensure the jar file exists
        if not os.path.isfile(self._jarfile):
            raise IOError('The required JAR file is not present. Reinstall jcompoundmapper-pywrapper.')
        # Define internal parameters
        self.fp_name = fingerprint if isinstance(fingerprint, str) else fingerprint.name
        self.fp_params = params

    def calculate(self, mols: Iterable[Chem.Mol], nbits: int = 1024, show_banner: bool = True, njobs: int = 1,
                  chunksize: Optional[int] = 1000) -> pd.DataFrame:
        """Calculate molecular fingerprints.

        :param mols: RDKit molecules for which fingerprints should be calculated
        :param nbits: size of fingerprints
        :param show_banner: If True, show notice on fingerprint usage
        :param njobs: number of concurrent processes
        :param chunksize: number of molecules to be processed by a process; ignored if njobs is 1
        :return: a pandas DataFrame containing the fingerprint values
        """
        if self.fp_name == 'MACCS':
            nbits = 166
        if show_banner:
            self._show_banner()
        # Parallelize should need be
        if njobs > 1:
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._multiproc_calculate, list(chunk), nbits)
                           for chunk in more_itertools.batched(mols, chunksize)
                           ]
            return pd.concat([future.result()
                              for future in futures]
                             ).reset_index(drop=True).fillna(0).astype(int)
        # Single process
        return self._calculate(list(mols), nbits=nbits)

    def _show_banner(self):
        """Print info message for citing."""
        print("""Fingerprinting algorithms for your use.
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

    def _prepare_command(self, mols: List[Chem.Mol], nbits: int = 1024) -> str:
        """Create the JCompoundMapper command to be run to obtain molecular fingerprints.

        :param mols: molecules to obtained molecular fingerprint of
        :param nbits: the size of the output fingerprints
        :return: The command to run.
        """
        # 1) Ensure JRE is accessible
        with self.lock:
            self._java_path = install_java()
        # 2) Create temp SD v2k file
        self._tmp_sd = mktempfile('molecules_v2k.sd')
        self._out = mktempfile('jcm_output.tsv')
        self._skipped = []
        try:
            block = BlockLogs()
            writer = Chem.SDWriter(self._tmp_sd)
            # Ensure V2000 as jCompoundMapper cannot properly process v3000
            writer.SetForceV3000(False)
            for i, mol in enumerate(mols):
                if mol is not None and isinstance(mol, Chem.Mol):
                    if mol.GetNumAtoms() > 999:
                        raise ValueError('Cannot calculate fingerprint for molecules with more than 999 atoms.')
                    # Does molecule lack hydrogen atoms?
                    if needsHs(mol):
                        warnings.warn('Molecule lacks hydrogen atoms: this might affect the value of calculated fingerprint')
                    # 3D fingerprint
                    if '3D' in self.fp_name:
                        confs = list(mol.GetConformers())
                        if not (len(confs) > 0 and confs[-1].Is3D()):
                            raise ValueError('Cannot calculate the 3D fingerprint of a conformer-less molecule')
                    writer.write(mol)
                else:
                    self._skipped.append(i)
            writer.close()
            del block
        except ValueError as e:
            # Free resources and raise error
            writer.close()
            del block
            os.remove(self._tmp_sd)
            raise e from None
        # 3) Create command
        java_path = install_java()
        if self.fp_params is not None:
            command_params = (f'-a {self.fp_params.atom_type.name} ' if self.fp_params.atom_type is not None else '') + \
                             (f'-d {self.fp_params.depth} ' if self.fp_params.depth is not None else '') + \
                             (f'-d {self.fp_params.dist_cutoff} ' if self.fp_params.dist_cutoff is not None else '') + \
                             (f'-k {str(self.fp_params.arom_flag).lower()} ' if self.fp_params.arom_flag is not None else '') + \
                             (f'-s {self.fp_params.stretch_factor} ' if self.fp_params.stretch_factor is not None else '')
        else:
            command_params = ''
        command = f"{java_path} -Djava.awt.headless=true -jar {self._jarfile} -f {self._tmp_sd} -ff LIBSVM_SPARSE " \
                  f"-hs {nbits} {command_params} -o {self._out}"
        return command

    def _cleanup(self, sd: bool = True, output: bool = True) -> None:
        """Cleanup resources used for calculation."""
        # Remove temporary files
        if sd:
            os.remove(self._tmp_sd)
        if output:
            os.remove(self._out)

    def _run_command(self, command: str, nbits: int = 1024) -> pd.DataFrame:
        """Run the jCompoundMapper command.

        :param command: The command to be run.
        :param nbits: size of output fingerprints
        """
        process = subprocess.run(command.split(), stdout=subprocess.DEVNULL)
        if process.returncode == 0:
            values = read_libsvmsparse(self._out, nbits)
        else:
            self._cleanup()
            raise ValueError('A conflicting set of parameters was used. Please fix or use the defaults.')
        values = pd.DataFrame(values, columns=[f'{self.fp_name}_{i}' for i in range(1, nbits + 1)])
        return values

    def _calculate(self, mols: List[Chem.Mol], nbits: int = 1024) -> pd.DataFrame:
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
                                              values=[np.NaN] * len(results.columns),
                                              axis=0),
                                    columns=results.columns)
                       )
        results = (results.apply(pd.to_numeric, errors='coerce', axis=1)
                          .fillna(0)
                          .convert_dtypes()
                   )
        return results

    def _multiproc_calculate(self, mols: List[Chem.Mol], nbits: int = 1024) -> pd.DataFrame:
        """Calculate jCompoundMapper fingerprints in thread-safe manner.

        :param mols: RDKit molecules for which jCompoundMapper fingerprints should be calculated
        :param nbits: size of the fingerprint
        :return: a pandas DataFrame containing the jCompoundMapper fingerprint values
        """
        # Copy self instance to make thread safe
        jcm = deepcopy(self)
        # Run copy
        result = jcm.calculate(mols, nbits=nbits, show_banner=False, njobs=1)
        return result
