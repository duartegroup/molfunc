from typing import Sequence, Optional
from molfunc.fragments import _frag_smiles_to_xyz_filenames
from molfunc_ext import (c_print_combined_from_names,
                         c_print_combined_from_xyz_filenames,
                         c_print_all_combined_molecules)


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
    _validate_xyz_filename(core_xyz_filename)

    if all(f is None for f in (frag_names, frag_smiles, frag_xyz_filenames)):
        raise ValueError('Combination requires at least some fragments!')

    if frag_smiles is not None:
        frag_xyz_filenames = _frag_smiles_to_xyz_filenames(frag_smiles)

    if frag_xyz_filenames and frag_names:
        raise ValueError('Cannot have fragment specified using structures and '
                         'names. Ordering may be ambiguous.')

    _decrement_atom_indices(atom_indices=atoms_to_del)

    # Cython passing strings python-> C requires encoding
    filename = bytes(f'{name}.xyz', 'utf-8')
    core_xyz_filename = bytes(core_xyz_filename, 'utf-8')

    if frag_names is not None:
        c_print_combined_from_names(filename,
                                    core_xyz_filename,
                                    atoms_to_del,
                                    [bytes(name, encoding='utf-8') for name in frag_names])
    else:
        c_print_combined_from_xyz_filenames(filename,
                                            core_xyz_filename,
                                            atoms_to_del,
                                            [bytes(fn, encoding='utf-8') for fn in frag_xyz_filenames])
    return None


def print_all_combined_molecules(core_xyz_filename: str,
                                 atoms_to_del:      Sequence[int],
                                 name:              str = 'combined_molecule'):
    """
    For a core given as a .xyz file generate all possible combinations of
    fragments defined in the library. Will generate N^m structures in the
    same .xyz file, where N is the number of fragments in the library and
    m is the number of atoms to delete

    ---------------------------------------------------------------------------
    Arguments:

        core_xyz_filename (str): .xyz filename of the core

        atoms_to_del (sequence(int)): Atoms to delete & exchange for fragments
                                      in the core. Indexed from 1

        name (str): Name of the molecule, thus generated .xyz file. e.g. ethane
                    -> ethane.xyz
    """
    _validate_xyz_filename(core_xyz_filename)
    _decrement_atom_indices(atom_indices=atoms_to_del)

    c_print_all_combined_molecules(bytes(f'{name}.xyz', 'utf-8'),
                                   bytes(core_xyz_filename, 'utf-8'),
                                   atoms_to_del)
    return None


def _validate_xyz_filename(xyz_filename: str):
    """Validate a .xyz filename"""

    if not xyz_filename.endswith('.xyz'):
        raise ValueError('Core .xyz filename must end with ".xyz". '
                         f'Had: {xyz_filename}')

    return None


def _decrement_atom_indices(atom_indices):

    for i, atom_idx in enumerate(atom_indices):
        if atom_idx <= 0:
            raise ValueError(f'Cannot have an atom index of {atom_idx}, '
                             f'must be indexed from 1')
        atom_indices[i] -= 1

    return None
