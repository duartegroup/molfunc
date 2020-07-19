from molfunc.molecules import FragmentMolecule
from molfunc.exceptions import MolFuncCritical
import os

here = os.path.dirname(os.path.abspath(__file__))
fragments_dir = os.path.join(here, 'fragments_lib')


class LibFragmentMolecule(FragmentMolecule):

    def __init__(self, name, filename):
        """Fragment with name aliases"""
        super().__init__(name=name, xyz_filename=filename)

        self.smiles = None
        self.aliases = []


def fragment_molecule_from_file(filename):
    """From a fragment in fragments_lib get a FragmentMolecule"""
    mol = LibFragmentMolecule(name=os.path.basename(filename).rstrip('.xyz'),
                              filename=filename)

    mol.smiles, mol.aliases = get_smiles_aliases(filename)

    return mol


def get_smiles_aliases(filename):
    """For a fragment molecule in a file get the aliases for it
    e.g. Me has aliased methyl and ch3"""

    for i, line in enumerate(open(filename, 'r')):
        # Aliases on the second line in the file
        if i == 1 and len(line.split()) == 2:
            # Line should be in the format: "SMILES alias1,alias2,alias3"
            smiles, aliases_string = line.split()
            return smiles, aliases_string.split(',')

    raise MolFuncCritical(f'Fragment molecule in {filename} had no aliases')


def get_fragment_molecule(name=None, smiles=None):
    """From a name e.g. Me get the corresponding FragmentMolecule"""
    if name is None and smiles is None:
        raise MolFuncCritical('Cannot get the fragment')

    # Iterate through all the fragments and return name or smiles matches
    for fragment_molecule in fragments:
        if name is not None and name.lower() in fragment_molecule.aliases:
            return fragment_molecule

        if smiles is not None and smiles == fragment_molecule.smiles:
            return fragment_molecule

    return None


# Populated when imported...
xyz_filepaths = [os.path.join(fragments_dir, fn)
                 for fn in os.listdir(fragments_dir) if fn.endswith('.xyz')]

# From all the xyz files populate fragments
fragments = [fragment_molecule_from_file(fn) for fn in xyz_filepaths]

# List of all aliases for fast checking if it exists
all_aliases = []
for fragment in fragments:
    all_aliases += fragment.aliases

all_smiles = [fragment.smiles for fragment in fragments]
