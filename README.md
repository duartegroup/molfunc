[![Build Status](https://travis-ci.org/duartegroup/molfunc.svg?branch=master)](https://travis-ci.org/duartegroup/molfunc) [![codecov](https://codecov.io/gh/duartegroup/molfunc/branch/master/graph/badge.svg)](https://codecov.io/gh/duartegroup/molfunc) [![PyPI version](https://badge.fury.io/py/molfunc.svg)](https://badge.fury.io/py/molfunc) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

![alt text](molfunc/common/example.png)

# molfunc

## About
**molfunc** is a Python tool for *func*tionalisation of 3D molecules. Given a
[.xyz](https://en.wikipedia.org/wiki/XYZ_file_format) file molecular functionalisation is performed
by specifying monovalent atoms in the structure to swap for a set of fragments given as their
corresponding SMILES strings. The energy of the combined molecule is minimised with purely rigid body rotations. 
Possible use cases include catalyst functionalisation, ligand modification and combinatorial molecule generation.

***
## Installation

If the requirements (rdkit, numpy, scipy, networkx) are already satisfied:
```
pip install molfunc
```

Otherwise, clone this repository and `cd` into the top level molfunc directory:
```
git clone https://github.com/duartegroup/molfunc.git
cd molfunc/
```
install the Python dependencies using conda  ([anaconda](https://www.anaconda.com/distribution/) or 
[miniconda](https://docs.conda.io/en/latest/miniconda.html)) 
and install:

```
conda config --append channels conda-forge
conda install --file requirements.txt
python setup.py install
```

***
## Usage
A minimal example to convert PH<sub>3</sub> to PMe<sub>3</sub>

```python
from molfunc import CoreMolecule, CombinedMolecule

ph3 = CoreMolecule(xyz_filename='examples/PH3.xyz', atoms_to_del=[2, 3, 4])
pme3 = CombinedMolecule(core_mol=ph3, frag_smiles='C[*]', name='PMe3')
pme3.print_xyz_file()
```

**molfunc** can also be used from the command line

```
molfunc examples/PH3.xyz -a 2 3 4 -s C[*]
```

where in both cases all hydrogen atoms are swapped for methyls. See *examples/* 
for more examples.
