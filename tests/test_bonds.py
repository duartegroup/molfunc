from molfunc.bonds import get_avg_bond_length


def test_bonds():

    # CH bond length ~ 1.0 Å
    assert 0.9 < get_avg_bond_length('C', 'H') < 1.1
    assert 0.9 < get_avg_bond_length('H', 'C') < 1.1

    # R-X bond length ~ 1.2 Å (because of RDKit)
    assert 1.1 < get_avg_bond_length('C', 'R') < 1.3

    # Unknown should be ~ 1.5
    assert 1.4 < get_avg_bond_length('C', 'X') < 1.6

    # Valid atoms but unknown interaction should be sensible
    assert 2.0 < get_avg_bond_length('Os', 'Sn') < 3.0
