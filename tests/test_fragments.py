from molfunc.fragments import fragments, all_aliases, all_smiles
from molfunc.fragments import LibFragmentMolecule
from molfunc.fragments import get_fragment_molecule
from molfunc.fragments import get_smiles_aliases
from molfunc.exceptions import MolFuncCritical
import pytest
import os

here = os.path.dirname(os.path.abspath(__file__))


def test_general_structure():

    assert len(fragments) > 0
    assert isinstance(fragments[0], LibFragmentMolecule)

    assert 'me' in all_aliases
    assert 'C[*]' in all_smiles


def test_get_fragment_molecule():

    methyl = get_fragment_molecule(name='me')
    assert methyl.n_atoms == 4

    methyl = get_fragment_molecule(smiles='C[*]')
    assert methyl.n_atoms == 4

    with pytest.raises(MolFuncCritical):
        # Cannot get a fragment with no name or SMILES string
        _ = get_fragment_molecule()

    not_a_valid_fragment = get_fragment_molecule(name='xxxx')
    assert not_a_valid_fragment is None


def test_fragment_from_file():

    benzene_file_path = os.path.join(here, 'data', 'benzene.xyz')

    with pytest.raises(MolFuncCritical):
        # This xyz file has no aliases, so should break making the fragment
        _ = get_smiles_aliases(filename=benzene_file_path)
