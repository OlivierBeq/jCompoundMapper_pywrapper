"""Tests for :mod:`jcompoundmapper_pywrapper.utils`."""

import os
import unittest

import numpy as np
from rdkit import Chem

from jcompoundmapper_pywrapper.utils import (
    get_java_in_dir,
    install_java,
    mktempdir,
    mktempfile,
    needsHs,
    read_libsvmsparse,
)


class TestNeedsHs(unittest.TestCase):
    """Tests for needsHs."""

    def test_needs_hs_without_explicit_hydrogens(self):
        mol = Chem.MolFromSmiles("CCO")
        self.assertTrue(needsHs(mol))

    def test_does_not_need_hs_with_explicit_hydrogens(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertFalse(needsHs(mol))


class TestTempFileHelpers(unittest.TestCase):
    """Tests for mktempdir/mktempfile."""

    def test_mktempdir_creates_writeable_dir(self):
        path = mktempdir()
        self.addCleanup(os.rmdir, path)
        self.assertTrue(os.path.isdir(path))

    def test_mktempfile_creates_writeable_file(self):
        path = mktempfile(suffix="foo.txt")
        self.addCleanup(os.remove, path)
        self.assertTrue(os.path.isfile(path))
        self.assertTrue(path.endswith("foo.txt"))


class TestReadLibsvmSparse(unittest.TestCase):
    """Tests for read_libsvmsparse."""

    def test_read_libsvmsparse(self):
        path = mktempfile(suffix=".tsv")
        self.addCleanup(os.remove, path)
        with open(path, "w") as fh:
            fh.write("mol1 1:1 3:1\n")
            fh.write("mol2 2:1\n")
        values = read_libsvmsparse(path, length=4)
        expected = np.array([[1, 0, 1, 0], [0, 1, 0, 0]], dtype=float)
        np.testing.assert_array_equal(values, expected)

    def test_read_libsvmsparse_empty_bits(self):
        path = mktempfile(suffix=".tsv")
        self.addCleanup(os.remove, path)
        with open(path, "w") as fh:
            fh.write("mol1\n")
        values = read_libsvmsparse(path, length=3)
        np.testing.assert_array_equal(values, np.zeros((1, 3)))


class TestJavaInstall(unittest.TestCase):
    """Tests for install_java/get_java_in_dir."""

    def test_install_java_returns_existing_path(self):
        path = install_java()
        self.assertIsNotNone(path)
        self.assertTrue(os.path.isfile(path))

    def test_get_java_in_dir_missing_version_returns_none(self):
        search_dir = mktempdir()
        self.addCleanup(os.rmdir, search_dir)
        path = get_java_in_dir(search_dir, version=999)
        self.assertIsNone(path)


if __name__ == "__main__":
    unittest.main()
