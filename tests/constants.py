"""Constants for unit tests."""

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

DATA_DIR = Path(__file__).parent / "data"


def load_diverse_subset(n: int = 250, with_3d: bool = False) -> list[Chem.Mol]:
    """Load a slice of the diverse 2000-molecule subset for broader test coverage.

    :param n: number of molecules to load from the front of the subset file
    :param with_3d: if True, embed a 3D conformer on every molecule (slower)
    :return: a list of RDKit molecules with explicit hydrogens
    """
    with open(DATA_DIR / "diverse_subset_2k.smi") as fh:
        smiles = [line.strip() for line in fh.readlines()[:n] if line.strip()]
    mols = [Chem.AddHs(mol) for smi in smiles if (mol := Chem.MolFromSmiles(smi)) is not None]
    if with_3d:
        for mol in mols:
            AllChem.EmbedMolecule(mol, randomSeed=42)
    return mols


MOLECULES = {
    "CHEMBL1560279":
        Chem.MolFromSmiles("CCN(CC)C(=O)[n+]1ccc(OC)cc1"),
    "erlotinib":
        Chem.MolFromSmiles("n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C"),
    "midecamycin":
        Chem.MolFromSmiles(
            "CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O"
            "[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)"
            "(C)O)N(C)C)O)CC=O)C)O)C"
        ),
    "selenofolate":
        Chem.MolFromSmiles("C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N"),
    "CHEMBL457061":
        Chem.MolFromSmiles("CCN(CC)C(=O)C1=CC=C(S1)C2=C3C=CC(=[N+](C)C)C=C3[Se]C4=C2C=CC(=C4)N(C)C"),
    "lomitapide":
        Chem.MolFromSmiles(
            "C1CN(CCC1NC(=O)C2=CC=CC=C2C3=CC=C(C=C3)C(F)(F)F)CCCCC4(C5=CC=CC=C5C6=CC=CC=C64)"
            "C(=O)NCC(F)(F)F"
        ),
    "cisplatin":
        Chem.MolFromSmiles("N.N.Cl[Pt]Cl"),
}

# Add hydrogens
for key, value in MOLECULES.items():
    MOLECULES[key] = Chem.AddHs(value)

# Obtain 3D conformer
for _, mol in MOLECULES.items():
    AllChem.EmbedMolecule(mol)

# Small slice of the diverse 2000-molecule subset, used for parallelism/broader coverage tests.
DIVERSE_SUBSET_SMALL = load_diverse_subset(40)
