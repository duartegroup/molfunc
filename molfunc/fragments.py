from typing import Sequence


def _frag_smiles_to_xyz_filenames(smiles_list: Sequence[str]) -> Sequence[str]:
    """
    Convert a list of SMILES strings of the form ['[*]C', '[*]Br'] to a
    list of .xyz filenames in the correct format

    Arguments:
        smiles_list (sequence(str)):

    Returns:
        (sequence(str)): xyz filenames

    Raises:
        ModuleNotFoundError: If autodE is not found

        ValueError: For invalid SMILES strings
    """
    try:
        from autode import Molecule

    except ModuleNotFoundError:
        raise ModuleNotFoundError('Cannot construct fragments from SMILES '
                                  'strings without an autodE install. Try:'
                                  '\n conda install autode -c conda-forge')

    xyz_filenames = [f'{idx}.xyz' for idx, _ in enumerate(smiles_list)]

    for smiles, filename in zip(smiles_list, xyz_filenames):

        if '[Li]' in smiles:
            raise NotImplementedError('Unfortunately SMILES fragments '
                                      'containing Li cannot be built')
        if '[*]' not in smiles:
            raise ValueError('Fragment SMILES must contain a [*]')

        mol = Molecule(smiles=smiles.replace('[*]', '[Li]'))
        for atom in mol.atoms:
            if atom.label == 'Li':
                atom.label = 'R'

        mol.print_xyz_file(filename=filename)

    return xyz_filenames

names = ["bz", "oph", "co2me", "bn", "oh", "ac", "tf", "ms", "me", "cn", "boc", "cl", "i", "h", "otbu", "ph", "ome", "tms", "et", "nme2", "nh2", "tbu", "co2et", "br", "ipr", "f", "no2", "ch2oh", "mes", "cf3"]
