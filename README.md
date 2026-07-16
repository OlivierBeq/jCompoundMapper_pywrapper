# 🧬 jCompoundMapper Python Wrapper

<!-- Badges -->
<div align="center">

[![PyPI version](https://img.shields.io/pypi/v/jcompoundmapper-pywrapper.svg)](https://pypi.org/project/jcompoundmapper-pywrapper/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/jcompoundmapper-pywrapper.svg)](https://pypi.org/project/jcompoundmapper-pywrapper/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/OlivierBeq/jCompoundMapper_pywrapper/actions/workflows/tests.yml/badge.svg)](https://github.com/OlivierBeq/jCompoundMapper_pywrapper/actions/workflows/tests.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
<br>
</div>

A simple and reliable Python wrapper for calculating [jCompoundMapper](https://doi.org/10.1186/1758-2946-3-3) molecular fingerprints. This library takes care of installing a matching Java runtime, dispatching molecules to the bundled jCMapperCLI executable, and collecting the results into a tidy pandas DataFrame — so you can stay in RDKit/pandas-land.

## ✨ Features

- 🧩 **17 fingerprint types** — DFS, ASP, AP2D, AT2D, AP3D, AT3D, CATS2D, CATS3D, PHAP2POINT2D, PHAP3POINT2D, PHAP2POINT3D, PHAP3POINT3D, ECFP, ECFPVariant, LSTAR, SHED, RAD2D, RAD3D and MACCS, straight from jCompoundMapper.
- ☕ **Zero Java setup** — automatically downloads, caches and reuses a matching JRE on first use; nothing to install by hand.
- ⚡ **Parallel by design** — spread the work across multiple CPU cores with configurable `njobs`/`chunksize`, each worker running its own single-core-pinned JVM.
- 🧯 **Never silently misaligned** — molecules that fail to parse are skipped explicitly and reinserted as `NaN` rows, never dropped without a trace.
- 📊 **pandas-native output** — results come back as a ready-to-use `DataFrame`, one row per molecule.
- 🎛️ **Configurable fingerprints** — override depth, distance cutoff, stretch factor, atom typing and aromaticity flags per fingerprint via `FpParamType`.

## ✍️ Copyright and Citation Notice

Olivier J. M. Béquignon is **neither** the copyright holder of jCompoundMapper **nor** responsible for it. The work carried out here concerns solely the Python wrapper.

### Citing

If you use this wrapper in your research, please cite the original jCompoundMapper publication in addition to this software package:

1. **Original jCompoundMapper paper:**
   > Hinselmann, G., Rosenbaum, L., Jahn, A., Fechner, N., Zell, A. (2011), jCompoundMapper: An open source Java library and command-line tool for chemical fingerprints. *Journal of Cheminformatics*, 3, 3.
   > [DOI: 10.1186/1758-2946-3-3](https://doi.org/10.1186/1758-2946-3-3)

2. **This wrapper:**
   > Béquignon, O. J. M. *jCompoundMapper_pywrapper: a Python wrapper for jCompoundMapper molecular fingerprints.* https://github.com/OlivierBeq/jCompoundMapper_pywrapper

## 📦 Installation

```bash
pip install jcompoundmapper-pywrapper
```

Or from source:

```bash
git clone https://github.com/OlivierBeq/jCompoundMapper_pywrapper.git
pip install ./jCompoundMapper_pywrapper
```

## 🛠️ Requirements

- Python 3.11+
- [RDKit](https://www.rdkit.org/docs/Install.html)

## 💡 Usage

### Getting started

```python
from jcompoundmapper_pywrapper import JCompoundMapper
from rdkit import Chem

smiles_list = [
    # erlotinib
    "n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C",
    # midecamycin
    "CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)(C)O)N(C)C)O)CC=O)C)O)C",
    # selenofolate
    "C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
    # cisplatin
    "N.N.Cl[Pt]Cl",
]
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

jcm = JCompoundMapper()
print(jcm.calculate(mols))
```

By default, the `DFS` fingerprint is calculated with 1024 bits. Any of the following fingerprints can be requested instead: `DFS`, `ASP`, `AP2D`, `AT2D`, `AP3D`, `AT3D`, `CATS2D`, `CATS3D`, `PHAP2POINT2D`, `PHAP3POINT2D`, `PHAP2POINT3D`, `PHAP3POINT3D`, `ECFP`, `ECFPVariant`, `LSTAR`, `SHED`, `RAD2D`, `RAD3D` or `MACCS` (always 166 bits, regardless of `nbits`).

```python
jcm = JCompoundMapper('ECFP')
print(jcm.calculate(mols, nbits=2048))

# or, using the Fingerprint enum

from jcompoundmapper_pywrapper import Fingerprint

jcm = JCompoundMapper(Fingerprint.ECFP)
print(jcm.calculate(mols, nbits=2048))
```

### 3D fingerprints

> ⚠️ Molecules with 3D conformers must be provided to calculate 3D fingerprints (i.e. `AP3D`, `AT3D`, `CATS3D`, `PHAP2POINT3D`, `PHAP3POINT3D` and `RAD3D`) — a `ValueError` is raised otherwise.
> ⚠️ A warning is raised if molecules lack explicit hydrogen atoms, since this may affect fingerprint values.

```python
from rdkit.Chem import AllChem

mols_3d = [Chem.AddHs(mol) for mol in mols]
_ = [AllChem.EmbedMolecule(mol) for mol in mols_3d]

jcm = JCompoundMapper('RAD3D')
print(jcm.calculate(mols_3d))
```

### Advanced fingerprint parameters

Each fingerprint ships with sensible defaults (search depth, distance cutoff, stretch factor, atom typing scheme and aromaticity flag). Custom parameters can be supplied via `FpParamType`:

```python
from jcompoundmapper_pywrapper import DEFAULT_FP_PARAMETERS

custom_ecfp_params = DEFAULT_FP_PARAMETERS['ECFP']
custom_ecfp_params.depth = 6

jcm = JCompoundMapper('ECFP', custom_ecfp_params)
print(jcm.calculate(mols))
```

### ⚡ Parallel processing

Speed things up by spreading molecules across several CPU cores. Each worker runs its own single-core-pinned JVM, so parallelism comes purely from the number of processes spawned — not from oversubscribing the host:

```python
jcm = JCompoundMapper('ECFP')
print(jcm.calculate(mols, njobs=8))
```

By default, molecules are auto-balanced evenly across `njobs` workers (`chunksize=None`), which minimizes JVM startup overhead while keeping every worker busy — the fastest setting for most workloads. A fixed `chunksize` can be provided instead if finer control is needed.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/OlivierBeq/jCompoundMapper_pywrapper/blob/master/LICENSE) file for details.

## 📚 API Documentation

```python
def calculate(mols, nbits=1024, show_banner=True, njobs=1, chunksize=None):
```

Calculates jCompoundMapper molecular fingerprints. Installs a matching JRE on first use if none is found.

#### Parameters

- ***mols  : Iterable[Chem.Mol]***
  RDKit molecule objects for which to obtain jCompoundMapper fingerprints.
- ***nbits  : int***
  Size of the fingerprints. Ignored for `MACCS`, which is always 166 bits.
- ***show_banner  : bool***
  Displays default notice about jCompoundMapper.
- ***njobs  : int***
  Number of concurrent processes used to calculate fingerprints in parallel; must not exceed the
  number of available CPU cores. Each spawned Java process is pinned to a single core
  (`-XX:ActiveProcessorCount=1`), since parallelism comes from spawning `njobs` OS processes
  rather than from letting each JVM oversubscribe the host's full core count.
- ***chunksize  : int | None***
  Number of molecules processed per worker process. If `None` (default), molecules are
  auto-balanced across all `njobs` workers so every worker gets work. Ignored if `njobs` is 1.
- ***return_type  : pd.DataFrame***
  Pandas DataFrame containing jCompoundMapper fingerprint values, one row per molecule.

________________

```python
JCompoundMapper(fingerprint='DFS', params=None)
```

Wrapper to obtain a single type of molecular fingerprint from jCompoundMapper.

#### Parameters

- ***fingerprint  : str | Fingerprint***
  Name of the fingerprint to calculate, either as a string or a `Fingerprint` enum member.
- ***params  : FpParamType | None***
  Custom fingerprint parameters. If `None` (default), the fingerprint's entry in
  `DEFAULT_FP_PARAMETERS` is used.

________________

```python
FpParamType(depth, dist_cutoff, stretch_factor, atom_type, arom_flag)
```

Parameter set of a jCompoundMapper fingerprint.

#### Attributes

- ***depth  : int | None***
  Search depth; only meaningful for 2D fingerprints (`DFS`, `ASP`, `AP2D`, `AT2D`, `CATS2D`,
  `PHAP2POINT2D`, `PHAP3POINT2D`, `ECFP`, `ECFPVariant`, `LSTAR`, `SHED`, `RAD2D`, `MACCS`).
- ***dist_cutoff  : float | None***
  Distance cutoff; only meaningful for 3D fingerprints (`AP3D`, `AT3D`, `CATS3D`,
  `PHAP2POINT3D`, `PHAP3POINT3D`, `RAD3D`).
- ***stretch_factor  : float | None***
  Stretch factor; only meaningful for 3D fingerprints.
- ***atom_type  : AtomType***
  Atom typing scheme used to seed the fingerprint.
- ***arom_flag  : bool***
  Whether aromaticity perception is used. Default: `False`.
