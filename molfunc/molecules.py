from typing import Sequence, Optional
from molfunc_ext import (print_combined_molecule_from_names,
                         print_combined_molecule_from_xyz_filenames)


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
        mol.print_xyz_file(filename=filename)

    return xyz_filenames


def print_combined_molecule(core_xyz_filename:  str,
                            atoms_to_del:       Sequence[int],
                            frag_smiles:        Optional[Sequence[str]] = None,
                            frag_names:         Optional[Sequence[str]] = None,
                            frag_xyz_filenames: Optional[Sequence[str]] = None,
                            name:               str = 'combined_molecule'):
    """
    Print an .xyz file of the combined molecule from a core defined in an .xyz
    file and a set of fragments.

    ---------------------------------------------------------------------------
    Arguments:

        core_xyz_filename (str): .xyz filename of the core

        atoms_to_del (sequence(int)): Atoms to delete & exchange for fragments
                                      in the core. Indexed from 1

        frag_smiles (sequence(str)): Names of the SMILES strings for the
                                     fragments (e.g. [*]C), where the position
                                     to functionalise is denoted with an *

        frag_names (sequence(str)): Names of the fragments e.g. Me for a methyl
                                    fragment. Supports aliases so any of
                                    {CH3, Me, methyl} afford the same behaviour

        frag_xyz_filenames (sequence(str)): Names of .xyz files that contain
                                            dummy atoms (R 'atomic' symbol)
                                            that are added to the core

        name (str): Name of the molecule, thus generated .xyz file. e.g. ethane
                    -> ethane.xyz
    """
    if not core_xyz_filename.endswith('.xyz'):
        raise ValueError('Core .xyz filename must end with ".xyz". '
                         f'Had: {core_xyz_filename}')

    if all(f is None for f in (frag_names, frag_smiles, frag_xyz_filenames)):
        raise ValueError('Combination requires at least some fragments!')

    if frag_smiles is not None:
        frag_xyz_filenames = _frag_smiles_to_xyz_filenames(frag_smiles)

    if frag_xyz_filenames and frag_names:
        raise ValueError('Cannot have fragment specified using structures and '
                         'names. Ordering may be ambiguous.')

    # Convert the atom indexes from 1-> indexing to 0-> indexing
    for i, atom_idx in enumerate(atoms_to_del):
        if atom_idx <= 0:
            raise ValueError(f'Cannot have an atom index of {atom_idx}, '
                             f'must be indexed from 1')
        atoms_to_del[i] -= 1

    if frag_names is not None:
        print_combined_molecule_from_names(f'{name}.xyz',
                                           core_xyz_filename,
                                           atoms_to_del,
                                           frag_names)
    else:
        print_combined_molecule_from_xyz_filenames(f'{name}.xyz',
                                                   core_xyz_filename,
                                                   atoms_to_del,
                                                   frag_xyz_filenames)
    return None
