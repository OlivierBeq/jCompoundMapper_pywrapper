"""Utility functions."""

import glob
import os
import sys
import tempfile

import numpy as np
from jdk import _JRE_DIR
from jdk import install as _jre_install
from rdkit import Chem


def install_java(version: int = 11) -> str | None:
    """Install a Java Runtime Environment."""
    path = get_java_in_dir(_JRE_DIR, version)
    if path is None:
        # Could not find JRE, install it
        _jre_install(version, jre=True)
        path = get_java_in_dir(_JRE_DIR, version=version)
    return path


def get_java_in_dir(dir: str, version: int) -> str | None:
    """Recursively search the directory to find a JRE."""
    paths = glob.glob(os.path.join(dir, '**', 'bin',
                                  'java.exe' if sys.platform == "win32" else 'java'
                                  ), recursive=True)
    path = [path for path in paths if f'jdk-{version}' in path]
    if len(path):
        return os.path.abspath(path[0])
    return None


def mktempdir(suffix: str | None = None) -> str:
    """Return the path to a writeable temporary directory."""
    return tempfile.mkdtemp(suffix=suffix)


def mktempfile(suffix: str | None = None) -> str:
    """Return the path to a writeable temporary file."""
    file = tempfile.mkstemp(suffix=suffix)
    os.close(file[0])
    return file[1]


def needsHs(mol: Chem.Mol) -> bool:
    """Return if the molecule lacks hydrogen atoms or not.

    :param mol: RDKit Molecule
    :return: True if the molecule lacks hydrogens.
    """
    for atom in mol.GetAtoms():
        nHNbrs = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                nHNbrs += 1
        noNeighbors = False
        if atom.GetTotalNumHs(noNeighbors) > nHNbrs:
            return True
    return False


def read_libsvmsparse(filename: str, length: int) -> np.ndarray:
    """Read the content of a LIBSVM_SPARSE file.

    :param filename: LIBSVM_SPARSE file
    :param length: size of the fingerprints contained in the file
    """
    # Read raw data
    with open(filename) as fh:  # file handle
        data = fh.readlines()
    results = []
    for result in data:
        # Get all bits set
        bits_set = result.split()[1:]
        # Create empty fp
        fp = np.zeros((length,))
        # Fill fp
        for x in bits_set:
            index, value = x.split(':')
            index = int(index) - 1
            fp[index] = int(value)
        results.append(fp)
    # Stack all fps
    return np.vstack(results)
