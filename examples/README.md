
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
Multiple atoms can be modified with the same fragment from the command line with 
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

### 3. Random substitution 
Random substitutions can be made, for example

```python
from molfunc import CoreMolecule, CombinedMolecule
import numpy as np

fragments = ['CN([*])C',
             'N[*]',
             'CCCN[*]',
             'O[*]',
             'CO[*]',
             'CCO[*]',
             'C[*]',
             'C[Si](C)([*])C',
             'F[*]',
             'Cl[*]',
             'Br[*]',
             'I[*]',
             '[*]C(OCC)=O',
             '[*]C(F)(F)F',
             'N#C[*]',
             'O=[N+]([*])[O-]']
             
benzene = CoreMolecule(xyz_filename='examples/benzene.xyz', atoms_to_del=[7, 9, 11])

# Select 3 random smiles strings from the list of possible fragments
combined = CombinedMolecule(core_mol=benzene,
                            frag_smiles_list=np.random.choice(fragments, size=3),
                            name=f'benzene_random_subst')
combined.print_xyz_file()
```

#### 4. Fragments from files
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


