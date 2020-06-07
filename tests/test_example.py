from molfunc import CoreMolecule, CombinedMolecule
from molfunc.atoms import xyz_file_to_atoms
import os
here = os.path.dirname(os.path.abspath(__file__))


def test_example():

    benzene_file_path = os.path.join(here, 'data', 'benzene.xyz')
    benzene = CoreMolecule(xyz_filename=benzene_file_path,
                           atoms_to_del=[7])

    toluene = CombinedMolecule(core_mol=benzene, frag_smiles='C[*]',
                               name='toluene')
    toluene.print_xyz_file()

    assert os.path.exists('toluene.xyz')
    gen_toluene_atoms = xyz_file_to_atoms('toluene.xyz')
    os.remove('toluene.xyz')

    # Toluene has 15 atoms
    assert len(gen_toluene_atoms) == 15
