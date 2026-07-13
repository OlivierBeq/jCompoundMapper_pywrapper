"""Tests for molecular fingerprints."""

import unittest
from unittest.mock import patch

from rdkit import Chem

from jcompoundmapper_pywrapper import DEFAULT_FP_PARAMETERS, FpParamType, JCompoundMapper
from tests.constants import MOLECULES


class TestFingerprints(unittest.TestCase):
    """Tests for jCompoundMapper_pywrapper molecular fingerprints."""
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


class TestJCompoundMapperValidation(unittest.TestCase):
    """Tests for input validation and error paths."""

    def test_invalid_fingerprint_name_raises(self):
        with self.assertRaises(ValueError):
            JCompoundMapper("NOT_A_FINGERPRINT")

    def test_invalid_params_type_raises(self):
        with self.assertRaises(ValueError):
            JCompoundMapper("DFS", params="not-a-param-object")

    def test_missing_jar_raises_ioerror(self):
        with patch.object(JCompoundMapper, "_jarfile", "/nonexistent/path/to.jar"):
            with self.assertRaises(IOError):
                JCompoundMapper("DFS")

    def test_atom_count_over_limit_raises(self):
        jcm = JCompoundMapper("DFS")
        big_mol = Chem.AddHs(Chem.MolFromSmiles("C" * 1000))
        with self.assertRaises(ValueError):
            jcm.calculate([big_mol], show_banner=False)

    def test_3d_fingerprint_without_conformer_raises(self):
        jcm = JCompoundMapper("AT3D")
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        with self.assertRaises(ValueError):
            jcm.calculate([mol], show_banner=False)

    def test_missing_hydrogens_warns(self):
        jcm = JCompoundMapper("DFS")
        mol = Chem.MolFromSmiles("CCO")
        with self.assertWarns(UserWarning):
            jcm.calculate([mol], show_banner=False)

    def test_skipped_molecule_yields_nan_free_row(self):
        jcm = JCompoundMapper("DFS")
        mols = [MOLECULES["erlotinib"], None, MOLECULES["cisplatin"]]
        values = jcm.calculate(mols, show_banner=False)
        self.assertEqual(values.shape[0], 3)
        self.assertFalse(values.isna().any().any())

    def test_custom_params_are_used(self):
        params = FpParamType(**vars(DEFAULT_FP_PARAMETERS["ECFP"]))
        params.depth = 4
        jcm = JCompoundMapper("ECFP", params)
        values = jcm.calculate([MOLECULES["erlotinib"]], show_banner=False)
        self.assertEqual(values.shape, (1, 1024))

    def test_subprocess_failure_surfaces_stderr(self):
        jcm = JCompoundMapper("DFS")
        command = jcm._prepare_command([MOLECULES["erlotinib"]])
        # Corrupt the jar argument so the JVM fails and stderr is captured.
        command[command.index("-jar") + 1] = "/nonexistent/bad.jar"
        with self.assertRaises(ValueError) as ctx:
            jcm._run_command(command, nbits=1024)
        self.assertTrue(str(ctx.exception))

    def test_show_banner_logs(self):
        jcm = JCompoundMapper("DFS")
        with self.assertLogs("jcompoundmapper_pywrapper.jcompoundmapper", level="INFO"):
            jcm.calculate([MOLECULES["erlotinib"]], show_banner=True)
