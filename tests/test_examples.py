from molfunc import CoreMolecule, CombinedMolecule
from scipy.spatial import distance_matrix
import numpy as np
import os

here = os.path.dirname(os.path.abspath(__file__))


def coordinates_are_resonable(coords):
    """Check that there are no very short or very long pairwise distances"""
    dist_mat = distance_matrix(coords, coords)
    return 0.8 < np.min(dist_mat + np.identity(len(coords))) < 5.0


def test_benzene_func():

    benzene_xyz_path = os.path.join(here, 'data', 'benzene.xyz')

    # Every hydrogen atom in benzene should be able to be functionalised
    for atom_to_del in [7, 8, 9, 10, 11, 12]:

        benzene = CoreMolecule(name='benzene',
                               xyz_filename=benzene_xyz_path,
                               atoms_to_del=[atom_to_del])

        toluene = CombinedMolecule(benzene, frag_smiles='C[*]')
        assert toluene.n_atoms == 15

    # Likewise if the atoms are reordered so the H atoms are first
    benzene_xyz_path = os.path.join(here, 'data', 'benzene_reordered.xyz')

    for atom_to_del in [1, 2, 3, 4, 5, 6]:

        benzene = CoreMolecule(name='benzene',
                               xyz_filename=benzene_xyz_path,
                               atoms_to_del=[atom_to_del])

        toluene = CombinedMolecule(benzene, frag_smiles='C[*]')
        assert toluene.n_atoms == 15


def test_dmhp():

    ph3 = CoreMolecule(xyz_filename=os.path.join(here, 'data', 'PH3.xyz'),
                       atoms_to_del=[2, 3, 4])

    dimphp = CombinedMolecule(core_mol=ph3,
                              frag_smiles_list=['C[*]', 'C[*]', '[*]C1=CC=CC=C1'],
                              name='dimphp')

    # Check that the geometry is fairly sensible
    assert coordinates_are_resonable(coords=dimphp.get_coordinates())

    dimphp.make_graph()
    assert dimphp.graph.number_of_edges() == 20


def test_pme3():

    ph3 = CoreMolecule(xyz_filename=os.path.join(here, 'data', 'PH3.xyz'),
                       atoms_to_del=[2, 3, 4])

    pme3 = CombinedMolecule(core_mol=ph3, frag_smiles='C[*]', name='PMe3')
    assert coordinates_are_resonable(coords=pme3.get_coordinates())

