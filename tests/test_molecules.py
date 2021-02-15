from molfunc.molecules import Molecule
from molfunc import FragmentMolecule, CombinedMolecule, CoreMolecule
from molfunc.atoms import Atom
from molfunc.exceptions import *
from scipy.spatial import distance_matrix
import numpy as np
import pytest
import os

here = os.path.dirname(os.path.abspath(__file__))


def coordinates_are_resonable(coords):
    """Check that there are no very short or very long pairwise distances"""
    dist_mat = distance_matrix(coords, coords)
    return 0.8 < np.min(dist_mat + np.identity(len(coords))) < 5.0


# Methane atoms
atoms = [Atom('C', 1.24788, 0.56457, -1.79703),
         Atom('H', 2.35728, 0.56457, -1.79703),
         Atom('H', 0.87808, 0.86789, -2.79804),
         Atom('H', 0.87807, 1.27982, -1.03385),
         Atom('H', 0.87807, -0.45398, -1.55920)]


def atoms_are_aliphatic_hydrocarbon(atoms_list):
    """Fully saturated hydrocarbon atoms"""

    for atom in atoms_list:
        if atom.label == 'H':
            assert atom.valence == 1       # H atoms are bonded to only C

        if atom.label == 'C':
            assert atom.valence == 4       # C atom is bonded to 4 Hs

    return None


def test_molecule_from_smiles():

    mol = Molecule(name='methane',
                        smiles='C')

    assert mol.n_atoms == 5
    assert mol.atoms is not None
    assert mol.graph is not None
    assert mol.name == 'methane'

    assert mol.graph.number_of_edges() == 4   # 4 CH bonds
    assert mol.graph.number_of_nodes() == 5   # 5 atoms

    atoms_are_aliphatic_hydrocarbon(mol.atoms)


def test_molecule_from_xyz_file():

    mol = Molecule(name='hydrocarbon',
                   xyz_filename=os.path.join(here, 'data', 'hydrocarbon.xyz'))

    assert mol.n_atoms == 20
    assert mol.atoms is not None
    assert mol.graph is not None
    assert mol.name == 'hydrocarbon'

    assert mol.graph.number_of_edges() == 19
    assert mol.graph.number_of_nodes() == 20

    atoms_are_aliphatic_hydrocarbon(mol.atoms)


def test_molecule_from_atoms():

    water_atoms = [Atom('O', -0.2594, -1.7353, -0.2594),
                   Atom('H',  0.3618, -1.3220,  0.3618),
                   Atom('H', -0.6950, -0.9849, -0.6950)]

    mol = Molecule(name='water', atoms=water_atoms)
    assert mol.n_atoms == 3
    assert mol.atoms is not None
    assert mol.graph is not None

    assert mol.atoms[0].valence == 2        # O atom is bonded to 2 hydrogens
    assert mol.atoms[1].valence == 1        # H atom is bonded to 1 oxygen
    assert mol.atoms[2].valence == 1        # H atom is bonded to 1 oxygen


def test_molecule_no_atoms():

    mol = Molecule()
    assert mol.n_atoms == 0
    assert mol.atoms is None
    assert mol.graph is None

    with pytest.raises(NoAtomsInMolecule):
        mol.print_xyz_file()

    with pytest.raises(NoAtomsInMolecule):
        mol.translate(vec=np.zeros(3))

    with pytest.raises(NoAtomsInMolecule):
        _ = mol.get_coordinates()


def test_molecule_bad_smiles():

    with pytest.raises(RDKitFailed):
        # Initialisation with an invalid SMILES string
        _ = Molecule(smiles='XX')


def test_core_molecule():

    core_mol = CoreMolecule(atoms=atoms)
    # Core molecule with no atoms to delete should retain all atoms
    assert core_mol.n_atoms == 5

    core_mol = CoreMolecule(atoms=atoms, atoms_to_del=[2])
    assert core_mol.n_atoms == 4

    # atoms_to_del list is indexed from 1, so atom 0 is not valid
    with pytest.raises(DatomsNotValid):
        _ = CoreMolecule(atoms=atoms, atoms_to_del=[0])

    with pytest.raises(DatomsNotValid):
        _ = CoreMolecule(atoms=atoms, atoms_to_del=[6])

    # Cannot use a carbon atom in ethane as the atom to delete
    xyz_path = os.path.join(here, 'data', 'ethane.xyz')

    with pytest.raises(DatomsNotValid):
        _ = CoreMolecule(name='ethane', xyz_filename=xyz_path, atoms_to_del=[1])

    # Cannot replace the same atom more than once
    with pytest.raises(DatomsNotValid):
        _ = CoreMolecule(name='methane', smiles='C',
                         atoms_to_del=[2, 2])


def test_bad_core_molecule():

    atoms_far_h = [Atom('C', 1.24788, 0.56457, -1.79703),
                   Atom('H', 9, 9, 9),
                   Atom('H', 0.87808, 0.86789, -2.79804),
                   Atom('H', 0.87807, 1.27982, -1.03385),
                   Atom('H', 0.87807, -0.45398, -1.55920)]

    # Core molecule with a atom to delete that doesn't have a nearest neighbour
    core_mol = CoreMolecule(atoms=atoms_far_h)

    with pytest.raises(DatomsNotValid):
        core_mol.get_datom_nn(datom_idx=1)


def test_fragment_molecule():

    mol = FragmentMolecule(smiles='C[*]')
    # Generating a fragment should delete the R atom
    assert mol.n_atoms == 4

    # Nearest neighbour to the R atom is C
    assert mol.nn_atom.label == 'C'

    ratoms = [Atom('C', 1.24788, 0.56457, -1.79703),
              Atom('R', 2.35728, 0.56457, -1.79703),
              Atom('H', 0.87808, 0.86789, -2.79804),
              Atom('H', 0.87807, 1.27982, -1.03385),
              Atom('H', 0.87807, -0.45398, -1.55920)]

    # Initialisation with atoms should also work
    mol = FragmentMolecule(atoms=ratoms)
    assert mol.n_atoms == 4
    assert mol.nn_atom.label == 'C'

    # Cannot generate a fragment from a structure where the R atom has > 1
    # bonds
    xyz_path = os.path.join(here, 'data', 'ethane_fragment.xyz')
    with pytest.raises(RAtomInvalidValence):
        _ = FragmentMolecule(xyz_filename=xyz_path)

    # Cannot generate a fragment from a structure with no R atoms
    xyz_path = os.path.join(here, 'data', 'ethane.xyz')
    with pytest.raises(RAtomNotFound):
        _ = FragmentMolecule(xyz_filename=xyz_path)


def test_fragment_molecule2():

    core_mol = CoreMolecule(atoms=atoms, atoms_to_del=[2])

    # Fragment molecules can also be initialised from [Fr] SMILES strings as
    # well as ['*']
    fragment_mol = FragmentMolecule(name='HOR', smiles='O[*]')
    assert fragment_mol.n_atoms == 2                            # OH atoms

    mol = CombinedMolecule(core_mol=core_mol, fragment=fragment_mol)
    assert mol.n_atoms == 6
    assert coordinates_are_resonable(coords=mol.get_coordinates())

    # Should also work directly from the SMILES string
    mol = CombinedMolecule(core_mol=core_mol, frag_smiles='O[*]')
    assert mol.n_atoms == 6
    assert coordinates_are_resonable(coords=mol.get_coordinates())


def test_bad_fragment_molecule():

    with pytest.raises(RAtomNotFound):
        _ = FragmentMolecule(atoms=[Atom('H')])


def test_fragment_molecule_closest():

    # R atom that is far away from the molecule is supported as the nn atom is
    # just the atom with the shortest distance
    xyz_path = os.path.join(here, 'data', 'ethane_fragment_far.xyz')
    mol = FragmentMolecule(xyz_filename=xyz_path)

    assert mol.nn_atom.label == 'C'


def test_combined_molecule():

    core_mol = CoreMolecule(atoms=atoms, atoms_to_del=[2])

    # Ethane
    mol = CombinedMolecule(core_mol=core_mol, frag_smiles='C[*]')
    assert coordinates_are_resonable(coords=mol.get_coordinates())
    assert mol.n_atoms == 8

    mol = CombinedMolecule(core_mol=core_mol, frag_smiles_list=['C[*]'])
    assert coordinates_are_resonable(coords=mol.get_coordinates())
    assert mol.n_atoms == 8

    fragment_mol = FragmentMolecule(smiles='C[*]')
    mol = CombinedMolecule(core_mol=core_mol, fragment=fragment_mol)
    assert coordinates_are_resonable(coords=mol.get_coordinates())
    assert mol.n_atoms == 8

    mol = CombinedMolecule(core_mol=core_mol, fragments=[fragment_mol])
    assert coordinates_are_resonable(coords=mol.get_coordinates())
    assert mol.n_atoms == 8

    # Delete two hydrogens from methae
    core_mol = CoreMolecule(atoms=atoms, atoms_to_del=[2, 3])

    # Propane
    mol = CombinedMolecule(core_mol=core_mol, frag_smiles='C[*]')
    assert coordinates_are_resonable(coords=mol.get_coordinates())
    assert mol.n_atoms == 11


def test_combined_molecule_no_core():

    core_mol = CoreMolecule(atoms=atoms)

    with pytest.raises(CombinationFailed):
        _ = CombinedMolecule(core_mol=core_mol, frag_smiles='C[*]')


def test_combined_molecule_diff_fragements():

    core_mol = CoreMolecule(atoms=atoms, atoms_to_del=[2])

    # Can't add two fragments to a molecule with one atom to delete
    with pytest.raises(CombinationFailed):
        _ = CombinedMolecule(core_mol=core_mol,
                             frag_smiles_list=['C[*]', 'C[*]'])

    # Similarly with two atoms to delete and a fragment smiles list of only
    # one fragment
    core_mol = CoreMolecule(atoms=atoms, atoms_to_del=[2, 3])

    with pytest.raises(CombinationFailed):
        _ = CombinedMolecule(core_mol=core_mol, frag_smiles_list=['C[*]'])


def test_combined_molecule_no_fragments():

    core_mol = CoreMolecule(atoms=atoms, atoms_to_del=[2])

    # Should be able to make a combined molecule from a core that has no
    # atoms to delete, but will add no atoms (Is this the clearest behaviour?)
    mol = CombinedMolecule(core_mol=core_mol)
    assert mol.n_atoms == 0


def test_combined_molecule_ortho_subst():

    methane = CoreMolecule(atoms=atoms, atoms_to_del=[2, 3])

    # Submitting two tBu groups onto methane should be possible (generate a
    # sensible structure if the optimisation includes the previously added
    # fragment
    subt = CombinedMolecule(core_mol=methane,
                            frag_smiles='CC(C)([*])C')

    assert coordinates_are_resonable(coords=subt.get_coordinates())


def test_combined_molecule_no_nn():

    methane = CoreMolecule(atoms=atoms)
    methane.datom_idxs = [2]

    with pytest.raises(CombinationFailed):
        _ = CombinedMolecule(core_mol=methane,
                             frag_smiles='C[*]')
