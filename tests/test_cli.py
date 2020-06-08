from subprocess import Popen
from molfunc.atoms import xyz_file_to_atoms
from molfunc.molecules import Molecule
from molfunc.molfunc import main
from scipy.spatial import distance_matrix
import numpy as np
import sys
import os

here = os.path.dirname(os.path.abspath(__file__))
xyz_path = os.path.join(here, 'data', 'benzene.xyz')


def test_cli():

    join = Popen(['molfunc', xyz_path, '-a', '7', '-s', 'C[*]'])
    join.wait()

    toluene_path = os.path.join(here, 'data', 'benzene_mod.xyz')
    assert os.path.exists(toluene_path)

    toluene_atoms = xyz_file_to_atoms(filename=toluene_path)
    assert len(toluene_atoms) == 15
    toluene = Molecule(atoms=toluene_atoms)

    # Geometry should be sensible.. i.e. min pairwise distance > 0.8 Å
    # and the maximum < 5 Å
    toluene_coords = toluene.get_coordinates()
    dist_mat = distance_matrix(toluene_coords, toluene_coords)

    assert 0.8 < np.min(dist_mat+np.identity(15)) < 5.0

    os.remove(toluene_path)


def test_molfunc():
    # For the example benzene -> toluene modification call the main molfunc
    sys.argv[1:] = [xyz_path, '-a', '7', '-s', 'C[*]']
    main()

    toluene_path = os.path.join(here, 'data', 'benzene_mod.xyz')
    assert os.path.exists(toluene_path)
    os.remove(toluene_path)
