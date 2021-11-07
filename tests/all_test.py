import os
import sys
import pytest
import numpy as np
from subprocess import Popen
from scipy.spatial import distance_matrix
from autode import Molecule
from molfunc.molfunc import main
from molfunc import (print_combined_molecule,
                     print_all_combined_molecules,
                     fragment_names)


here = os.path.dirname(os.path.abspath(__file__))
xyz_path = os.path.join(here, 'data', 'benzene.xyz')


def _xyz_file_is_reasonable(filename):
    """Is the structure within an .xyz file reasonable?"""

    if not os.path.exists(filename):
        return False  # Definitely doesn't without an existing file!

    mol = Molecule(filename)

    dist_mat = distance_matrix(mol.coordinates, mol.coordinates)
    dist_mat += np.identity(n=mol.n_atoms)   # Remove zero diagonal

    # No very short or very long distances
    return 0.7 < np.min(dist_mat) < 20.0


def test_cli():
    """Test command line interface"""

    join = Popen(['molfunc', xyz_path, '-a', '7', '-f', 'Me'])
    join.wait()

    toluene_path = os.path.join(here, 'data', 'benzene_mod.xyz')
    assert os.path.exists(toluene_path)

    toluene = Molecule(toluene_path)
    assert toluene.n_atoms == 15
    assert _xyz_file_is_reasonable(toluene_path)

    os.remove(toluene_path)


def test_main():
    # For the example benzene -> toluene modification call the main molfunc
    sys.argv[1:] = [xyz_path, '-a', '7', '-f', 'Me']
    main()

    toluene_path = os.path.join(here, 'data', 'benzene_mod.xyz')
    assert os.path.exists(toluene_path)
    os.remove(toluene_path)


def test_main_all():
    # For the example benzene -> toluene modification call the main molfunc
    sys.argv[1:] = [xyz_path, '-a', '7', '--all']
    main()

    toluene_path = os.path.join(here, 'data', 'benzene_mod.xyz')
    assert os.path.exists(toluene_path)
    os.remove(toluene_path)


def test_examples_fragment():

    print_combined_molecule(core_xyz_filename=xyz_path,
                            atoms_to_del=[7],
                            frag_names=['Me'],
                            name='toluene')

    assert _xyz_file_is_reasonable('toluene.xyz')
    os.remove('toluene.xyz')

    print_combined_molecule(core_xyz_filename=xyz_path,
                            atoms_to_del=[7, 9, 11],
                            frag_names=fragment_names[:3],
                            name='benzene_random_subst')

    assert _xyz_file_is_reasonable('benzene_random_subst.xyz')
    os.remove('benzene_random_subst.xyz')

    print_combined_molecule(core_xyz_filename=xyz_path,
                            atoms_to_del=[7],
                            frag_xyz_filenames=[os.path.join(here, 'data', 'methyl.xyz')],
                            name='toluene')
    assert _xyz_file_is_reasonable('toluene.xyz')
    os.remove('toluene.xyz')


def test_ok_smiles_fragment():

    print_combined_molecule(core_xyz_filename=xyz_path,
                            atoms_to_del=[7],
                            frag_smiles=['C[*]'],
                            name='toluene')

    assert _xyz_file_is_reasonable('toluene.xyz')
    os.remove('toluene.xyz')


def test_not_ok_smiles_fragment():

    with pytest.raises(Exception):
        # Fragment needs a [*] specification
        print_combined_molecule(core_xyz_filename=xyz_path,
                                atoms_to_del=[7],
                                frag_smiles=['C'])

    with pytest.raises(Exception):
        # Cannot have a Li atom in...
        print_combined_molecule(core_xyz_filename=xyz_path,
                                atoms_to_del=[7],
                                frag_smiles=['[Li][*]'])


def test_not_ok_core():

    with pytest.raises(Exception):
        # Must end in .xyz
        print_combined_molecule(core_xyz_filename=os.path.join(here, 'data', 'benzene'),
                                atoms_to_del=[7],
                                frag_names=['Me'])

    with pytest.raises(Exception):
        # Combination requires at least one fragment specification
        print_combined_molecule(core_xyz_filename=xyz_path,
                                atoms_to_del=[7])

    with pytest.raises(Exception):
        # Cannot have both a name and a .xyz files
        print_combined_molecule(core_xyz_filename=xyz_path,
                                atoms_to_del=[7],
                                frag_names=["Me"],
                                frag_xyz_filenames=[os.path.join(here, 'data', 'methyl.xyz')])

    with pytest.raises(Exception):
        # Indexing of atoms to user-facing code is from 1
        print_combined_molecule(core_xyz_filename=os.path.join(here, 'data', 'benzene.xyz'),
                                atoms_to_del=[0],
                                frag_names=['Me'])


def test_all_combination():
    """Test all possible combinations can be generated reasonably"""

    ph3_filepath = os.path.join(here, 'data', 'PH3.xyz')
    print_all_combined_molecules(ph3_filepath,
                                 atoms_to_del=[2],
                                 name='tmp')

    assert os.path.exists('tmp.xyz')

    xyz_lines = []
    for line in open('tmp.xyz', 'r'):

        xyz_lines.append(line)

        if len(line.split()) != 1 or len(xyz_lines) <= 1:
            continue

        with open('single_tmp.xyz', 'w') as single_xyz_file:
            for _line in xyz_lines[:-1]:
                print(_line, file=single_xyz_file, end='')

        assert _xyz_file_is_reasonable('single_tmp.xyz')
        xyz_lines = xyz_lines[-1:]

    os.remove('tmp.xyz')
    os.remove('single_tmp.xyz')
