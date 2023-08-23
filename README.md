[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Python wrapper for jCompoundMapper molecular fingerprints

Python wrapper to ease the calculation of jCompoundMApper molecular fingerprints.

## Installation

From source:

    git clone https://github.com/OlivierBeq/jcompoundmapper_pywrapper.git
    pip install ./jcompoundmapper_pywrapper

with pip:

```bash
pip install jcompoundmapper-pywrapper
```

### Get started

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
  "N.N.Cl[Pt]Cl"
]
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

jcm = JCompoundMapper()
print(jcm.calculate(mols))
```

One can use any of the following fingerprints: DFS, ASP, AP2D, AT2D, AP3D, AT3D, CATS2D, CATS3D, PHAP2POINT2D, PHAP3POINT2D, PHAP2POINT3D, PHAP3POINT3D,
ECFP, ECFPVariant, LSTAR, SHED, RAD2D, RAD3D or MACCS.

By default, 1024 bits fingerprints are generated, except for the 166 bits MACCS.

```python
jcm = JCompoundMapper('ECFP')
print(jcm.calculate(mols, 2048))

# or

from jCompoundMapper_pywrapper import Fingerprint

jcm = JCompoundMapper(Fingerprint.DFS)
print(jcm.calculate(mols, 2048))
```

:warning: Molecules with 3D conformers must be provided to calculate 3D fingerprints (i.e. CATS3D, PHAP2POINT3D, PHAP3POINT3D and RAD3D).

One can specify advanced parameters to the fingerprinter:

```python
from jCompoundMapper_pywrapper import DEFAULT_FP_PARAMETERS

custom_ecfp_params = DEFAULT_FP_PARAMETERS['ECFP']
custom_ecfp_params.depth = 6

jcm = JCompoundMapper('ECFP', custom_ecfp_params)
print(jcm.calculate(mols))
```
## Documentation

```python
def calculate(mols, nbits=1024, show_banner=True, njobs=1, chunksize=1000):
```

Default method to calculate jCompoundMapper fingerprints.

Parameters:

- ***mols  : Iterable[Chem.Mol]***  
  RDKit molecule objects for which to obtain jCompoundMapper fingerprints.
- ***nbits  : Union[int, List[int]]***  
  Size of the fingerprints.
- ***show_banner  : bool***  
  Displays default notice about jCompoundMapper.
- ***njobs  : int***  
  Maximum number of simultaneous processes.
- ***chunksize  : int***  
  Maximum number of molecules each process is in charge of.
