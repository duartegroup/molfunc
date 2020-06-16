from molfunc.atoms import xyz_file_to_atoms
from molfunc.atoms import smiles_to_atoms
from molfunc.atoms import Atom
from molfunc.exceptions import *
import pytest
import os
here = os.path.dirname(os.path.abspath(__file__))


def test_okay_xyz_file():

    atoms = xyz_file_to_atoms(os.path.join(here, 'data', 'hydrocarbon.xyz'))
    assert len(atoms) == 20

    for atom in atoms:
        assert isinstance(atom, Atom)

    # String representation should be concise
    assert len(str(atoms[0])) < 50


def test_broken_xyz_file():

    with pytest.raises(XYZfileDidNotExist):
        _ = xyz_file_to_atoms(filename='XXXXXX')

    with pytest.raises(NotXYZfile):
        path = os.path.join(here, 'data', 'wrong_extension_xyz_file.mol')
        _ = xyz_file_to_atoms(filename=path)

    with pytest.raises(XYZfileMalformatted):
        path = os.path.join(here, 'data', 'wrong_format_xyz_file.xyz')
        _ = xyz_file_to_atoms(filename=path)

    with pytest.raises(XYZfileMalformatted):
        path = os.path.join(here, 'data', 'no_data_xyz_file.xyz')
        _ = xyz_file_to_atoms(filename=path)


def test_smiles_generation():

    # Methane smiles should be 5 atoms
    assert len(smiles_to_atoms(smiles='C', n_confs=1)) == 5

    # Requesting more than 1 conformer of propane should generate a list of
    # lists
    conformer_list = smiles_to_atoms(smiles='CCC', n_confs=2)
    assert 1 <= len(conformer_list) <= 2

    # Propane should have 11 atoms in each conformer
    for conformer_atoms in conformer_list:
        assert len(conformer_atoms) == 11
