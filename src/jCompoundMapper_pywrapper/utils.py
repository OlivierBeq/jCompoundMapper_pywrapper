# -*- coding: utf-8

"""Utility functions."""

import sys
import os
import glob
import tempfile
import shutil
from pathlib import Path
from typing import Tuple

import numpy as np
from rdkit import Chem
from jdk import install as _jre_install, _JRE_DIR


def install_jre(version: int = 11):
    """Install a Java Runtime Environment."""
    path = get_jre_in_dir(_JRE_DIR)
    if len(path) == 0:
        # Could not find JRE, install it
        _ = _jre_install(version, jre=True)
        path = get_jre_in_dir(_JRE_DIR)
    return path


def make_temp_jre() -> Tuple[str, str]:
    """Copy the installed JRE to a temporary file.

    Allows multiprocessing capacities.

    :return: a tuple of (path to temp dir, path to the JRE)
    """
    # Path of JRE installation
    install_path = [path for path in Path(install_jre()).parents
                    if path.as_posix().endswith('jre') and not path.as_posix().endswith('.jre')][0]
    # Path of temp dir
    outdir = mktempdir()
    # Make copy into dir
    shutil.copytree(install_path, os.path.join(outdir, 'jre'))
    # Find JRE in temp dir
    path = get_jre_in_dir(outdir)
    return outdir, path


def get_jre_in_dir(dir: str):
    """Recursively search the directory to find a JRE."""
    path = glob.glob(os.path.join(dir, '**', 'server',
                                  'jvm.dll' if sys.platform == "win32" else 'libjvm.so'
                                  ), recursive=True)
    if len(path):
        return path[0]
    return None


def install_java(version: int = 11):
    """Install a Java Runtime Environment."""
    path = get_java_in_dir(_JRE_DIR)
    if path is None:
        # Could not find JRE, install it
        _ = _jre_install(version, jre=True)
        path = get_java_in_dir(_JRE_DIR)
    return path

def get_java_in_dir(dir: str):
    """Recursively search the directory to find a JRE."""
    path = glob.glob(os.path.join(dir, '**', 'bin',
                                  'java.exe' if sys.platform == "win32" else 'java'
                                  ), recursive=True)
    if len(path):
        return os.path.abspath(path[0])
    return None


def mktempdir(suffix: str = None) -> str:
    """Return the path to a writeable temporary directory."""
    dir = tempfile.mkdtemp(suffix=suffix)
    return dir


def mktempfile(suffix: str = None) -> str:
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
