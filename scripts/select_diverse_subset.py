"""One-off script to build a reproducible diverse molecule subset for tests.

Not shipped as part of the package. Samples a tractable candidate pool from a
large .smi file, then uses RDKit's MaxMin picker over Morgan fingerprints to
select a diverse subset. Run once; the output is committed under tests/data/.

Usage:
    python scripts/select_diverse_subset.py <input.smi> <output.smi> [--pool-size N] [--n-picks N] [--seed N]
"""

import argparse
import random

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

MAX_ATOMS = 999


def sample_candidate_pool(input_path: str, pool_size: int, seed: int) -> list[str]:
    with open(input_path) as fh:
        lines = [line.strip() for line in fh if line.strip()]
    random.seed(seed)
    return random.sample(lines, min(pool_size, len(lines)))


def filter_and_dedupe(smiles_list: list[str]) -> list[str]:
    RDLogger.DisableLog("rdApp.*")
    seen = set()
    kept = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None or mol.GetNumAtoms() == 0 or mol.GetNumAtoms() > MAX_ATOMS:
            continue
        canonical = Chem.MolToSmiles(mol)
        if canonical in seen:
            continue
        seen.add(canonical)
        kept.append(canonical)
    return kept


def pick_diverse(smiles_list: list[str], n_picks: int, seed: int) -> list[str]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [generator.GetFingerprint(Chem.MolFromSmiles(smi)) for smi in smiles_list]
    picker = MaxMinPicker()
    indices = picker.LazyBitVectorPick(fps, len(fps), min(n_picks, len(fps)), seed=seed)
    return [smiles_list[i] for i in indices]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_smi")
    parser.add_argument("output_smi")
    parser.add_argument("--pool-size", type=int, default=20_000)
    parser.add_argument("--n-picks", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    candidates = sample_candidate_pool(args.input_smi, args.pool_size, args.seed)
    valid = filter_and_dedupe(candidates)
    diverse = pick_diverse(valid, args.n_picks, args.seed)

    with open(args.output_smi, "w") as fh:
        fh.write("\n".join(diverse) + "\n")

    print(f"Sampled {len(candidates)} candidates, {len(valid)} valid/unique, picked {len(diverse)} diverse.")


if __name__ == "__main__":
    main()
