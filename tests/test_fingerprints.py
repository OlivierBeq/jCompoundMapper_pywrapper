# -*- coding: utf-8 -*-
"""Tests for molecular fingerprints."""

import unittest

from jcompoundmapper_pywrapper import JCompoundMapper, DEFAULT_FP_PARAMETERS
from tests.constants import MOLECULES


class TestFingerprints(unittest.TestCase):
    """Tests for CDK_pywrapper molecular fingerprints."""
    def setUp(self) -> None:
        """Load molecules."""
        self.molecules = list(MOLECULES.values())
        self.fingerprints = DEFAULT_FP_PARAMETERS.keys()

    def test_fingerprint(self):
        """Test the dimensions of the output fingerprint dataframe."""
        for fp_type in self.fingerprints:
            jcm = JCompoundMapper(fp_type)
            values = jcm.calculate(self.molecules, show_banner=False)
            if fp_type != "MACCS":
                self.assertEqual(values.shape, (len(MOLECULES), 1024))
                self.assertEqual(len(values.columns.unique().tolist()), 1024)
            else:
                self.assertTrue(values.shape == (len(MOLECULES), 166))
                self.assertEqual(len(values.columns.unique().tolist()), 166)
            self.assertFalse(values.isna().any().any())

    def test_fingerprint_2048bits(self):
        """Test the dimensions of the output fingerprint dataframe."""
        for fp_type in self.fingerprints:
            jcm = JCompoundMapper(fp_type)
            values = jcm.calculate(self.molecules, nbits=2048, show_banner=False)
            if fp_type != "MACCS":
                self.assertEqual(values.shape, (len(MOLECULES), 2048))
                self.assertEqual(len(values.columns.unique().tolist()), 2048)
            else:
                self.assertTrue(values.shape == (len(MOLECULES), 166))
                self.assertEqual(len(values.columns.unique().tolist()), 166)
            self.assertFalse(values.isna().any().any())

    def test_fingerprint_multithread(self):
        """Test the dimensions of the output fingerprint dataframes calculated by different processes."""
        for fp_type in self.fingerprints:
            jcm = JCompoundMapper(fp_type)
            values = jcm.calculate(self.molecules, show_banner=False,
                                   njobs=1, chunksize=1)
            if fp_type != "MACCS":
                self.assertEqual(values.shape, (len(MOLECULES), 1024))
                self.assertEqual(len(values.columns.unique().tolist()), 1024)
            else:
                self.assertTrue(values.shape == (len(MOLECULES), 166))
                self.assertEqual(len(values.columns.unique().tolist()), 166)
            self.assertFalse(values.isna().any().any())

    def test_fingerprint_2048bits_multithread(self):
        """Test the dimensions of the output fingerprint dataframe."""
        for fp_type in self.fingerprints:
            jcm = JCompoundMapper(fp_type)
            values = jcm.calculate(self.molecules, nbits=2048, show_banner=False,
                                   njobs=1, chunksize=1)
            if fp_type != "MACCS":
                self.assertEqual(values.shape, (len(MOLECULES), 2048))
                self.assertEqual(len(values.columns.unique().tolist()), 2048)
            else:
                self.assertTrue(values.shape == (len(MOLECULES), 166))
                self.assertEqual(len(values.columns.unique().tolist()), 166)
            self.assertFalse(values.isna().any().any())
