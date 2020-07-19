
### 0. Toluene 
To functionalise benzene with a methyl group generating toluene the SMILES 
string of the methyl group can be obtained from Chemdraw
 
![alt text](../molfunc/common/smiles_example.png)

To modify the 3D benzene structure in _examples/_ by replacing the first hydrogen atom (number 7) with a methyl
```
molfunc examples/benzene.xyz -a 7 -s C[*]
```
which will generate a benzene_mod.xyz file in the _examples/_ directory. Alternatively, in a Python script

```python
from molfunc import CoreMolecule, CombinedMolecule

benzene = CoreMolecule(xyz_filename='examples/benzene.xyz', atoms_to_del=[7])
toluene = CombinedMolecule(core_mol=benzene, frag_smiles='C[*]', name='toluene')
toluene.print_xyz_file()
```
which will generate toluene.xyz file in the current working directory.

![alt text](../molfunc/common/benzene_func.png)

### 1. _p_-Xylene
Multiple atoms can be modified with the same fragment, from the command line 
```
molfunc examples/benzene.xyz -a 7 10 -s C[*]
```
which swaps two para hydrogens to methyls and generates benzene_mod.xyz as _p_-xylene.


### 2. Dimethylphenylphosphine

Multiple modifications are available using the Python module. To generate dimethylphenylphosphine
from PH<sub>3</sub> 

```python
from molfunc import CoreMolecule, CombinedMolecule

ph3 = CoreMolecule(xyz_filename='examples/PH3.xyz', atoms_to_del=[2, 3, 4])
dimphp = CombinedMolecule(core_mol=ph3, 
                          frag_smiles_list=['C[*]', 'C[*]', '[*]C1=CC=CC=C1'], 
                          name='dimphp')
dimphp.print_xyz_file()
```

![alt text](../molfunc/common/ph3_func.png)

Note that the functionalisations are made in order so atoms 2 and 3 will be replaced with methyls 
and atom 4 with a phenyl.

### 3. Triphenylphoshpine

Common fragments are available in the fragment library and don't require an
RDKit install (and similarly in 5.)

```python
from molfunc import CoreMolecule, CombinedMolecule
from molfunc.fragments import get_fragment_molecule

ph3 = CoreMolecule(xyz_filename='examples/PH3.xyz', atoms_to_del=[2, 3, 4])
# Fragments are added with aliases so any of {Ph, phenyl, C6H5} are possible
ph = get_fragment_molecule(name='phenyl')

# Swap each hydrogen for the phenyl fragment
tpp = CombinedMolecule(core_mol=ph3, fragment=ph, name='TPP')
tpp.print_xyz_file()
```

### 4. Random substitution 
Random substitutions can be made, for example

```python
from molfunc import CoreMolecule, CombinedMolecule
from molfunc.fragments import fragments
import numpy as np

benzene = CoreMolecule(xyz_filename='examples/benzene.xyz', atoms_to_del=[7, 9, 11])

# Select 3 fragments from those in molfunc/fragments_lib/
combined = CombinedMolecule(core_mol=benzene,
                            fragments=np.random.choice(fragments, size=3),
                            name='benzene_random_subst')
combined.print_xyz_file()
```

#### 5. Fragments from files
Fragments can also be generated from xyz files. They must contain an 
'R' atom which will be swapped for the core molecule. For example, to
generate toluene 

```python
from molfunc import CoreMolecule, CombinedMolecule, FragmentMolecule

benzene = CoreMolecule(xyz_filename='examples/benzene.xyz', atoms_to_del=[7])
methyl = FragmentMolecule(name='methyl', xyz_filename='examples/methyl.xyz')

toluene = CombinedMolecule(core_mol=benzene, fragment=methyl, name='toluene')
toluene.print_xyz_file()
```


