"""Tests for parallel dispatch: chunk sizing, core-count validation, JVM pinning."""

import unittest
from unittest.mock import patch

from jcompoundmapper_pywrapper import JCompoundMapper
from jcompoundmapper_pywrapper.jcompoundmapper import _make_chunks
from tests.constants import DIVERSE_SUBSET_SMALL


class TestMakeChunks(unittest.TestCase):
    """Unit tests for the pure chunk-splitting helper (no Java involved)."""

    def test_no_molecules(self):
        self.assertEqual(_make_chunks([], njobs=4, chunksize=None), [])

    def test_auto_split_uses_all_workers_when_enough_molecules(self):
        mols = list(range(7))
        chunks = _make_chunks(mols, njobs=4, chunksize=None)
        self.assertEqual(len(chunks), 4)
        self.assertEqual(sorted(len(c) for c in chunks), [1, 2, 2, 2])
        self.assertEqual(sum(chunks, []), mols)

    def test_auto_split_fewer_molecules_than_workers(self):
        mols = list(range(3))
        chunks = _make_chunks(mols, njobs=8, chunksize=None)
        # Cannot use more workers than there are molecules: exactly len(mols) chunks.
        self.assertEqual(len(chunks), 3)
        self.assertEqual(sorted(len(c) for c in chunks), [1, 1, 1])

    def test_explicit_chunksize_is_authoritative(self):
        mols = list(range(7))
        chunks = _make_chunks(mols, njobs=4, chunksize=1000)
        # A too-large explicit chunksize still yields a single chunk (old batching behavior).
        self.assertEqual(len(chunks), 1)
        self.assertEqual(chunks[0], mols)

    def test_explicit_small_chunksize(self):
        mols = list(range(7))
        chunks = _make_chunks(mols, njobs=4, chunksize=3)
        self.assertEqual([len(c) for c in chunks], [3, 3, 1])


class TestNjobsValidation(unittest.TestCase):
    """Tests for njobs bounds checking, without invoking Java."""

    def test_njobs_zero_raises(self):
        jcm = JCompoundMapper("DFS")
        with self.assertRaises(ValueError):
            jcm.calculate(DIVERSE_SUBSET_SMALL[:2], njobs=0, show_banner=False)

    def test_njobs_exceeds_cpu_count_raises(self):
        jcm = JCompoundMapper("DFS")
        with patch("os.cpu_count", return_value=2):
            with self.assertRaises(ValueError):
                jcm.calculate(DIVERSE_SUBSET_SMALL[:2], njobs=3, show_banner=False)


class TestJvmCommandPinning(unittest.TestCase):
    """Tests that every Java invocation is pinned to a single core."""

    def test_prepare_command_pins_processor_count(self):
        jcm = JCompoundMapper("DFS")
        command = jcm._prepare_command(DIVERSE_SUBSET_SMALL[:2], nbits=64)
        self.addCleanup(jcm._cleanup)
        self.assertIsInstance(command, list)
        self.assertIn("-XX:ActiveProcessorCount=1", command)


class TestParallelMatchesSingleProcess(unittest.TestCase):
    """Ensures the njobs>1 code path is exercised and produces the same results as njobs=1."""

    def test_parallel_matches_single_process(self):
        mols = DIVERSE_SUBSET_SMALL[:12]
        single = JCompoundMapper("ECFP").calculate(mols, nbits=128, show_banner=False, njobs=1)
        parallel = JCompoundMapper("ECFP").calculate(
            mols, nbits=128, show_banner=False, njobs=4, chunksize=None
        )
        self.assertEqual(single.shape, parallel.shape)
        self.assertTrue((single.values == parallel.values).all())
        self.assertEqual(list(single.dtypes), list(parallel.dtypes))

    def test_parallel_uses_all_requested_workers(self):
        # With 12 molecules and njobs=4, _make_chunks must produce 4 chunks: verified directly
        # (avoids depending on OS-level process introspection, which is flaky in CI).
        mols = DIVERSE_SUBSET_SMALL[:12]
        chunks = _make_chunks(mols, njobs=4, chunksize=None)
        self.assertEqual(len(chunks), 4)


if __name__ == "__main__":
    unittest.main()
