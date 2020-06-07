from subprocess import Popen
from molfunc.atoms import xyz_file_to_atoms
from molfunc.molfunc import main
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
    os.remove(toluene_path)


def test_molfunc():
    # For the example benzene -> toluene modification call the main molfunc
    sys.argv[1:] = [xyz_path, '-a', '7', '-s', 'C[*]']
    main()

    toluene_path = os.path.join(here, 'data', 'benzene_mod.xyz')
    assert os.path.exists(toluene_path)
    os.remove(toluene_path)
