from typing import Sequence, Optional
from molfunc.fragments import _frag_smiles_to_xyz_filenames
from molfunc_ext import (print_combined_from_names,
                         print_combined_from_xyz_filenames)


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

    # Cython passing strings python-> C requires encoding
    filename = bytes(f'{name}.xyz', 'utf-8')
    core_xyz_filename = bytes(core_xyz_filename, 'utf-8')

    if frag_names is not None:
        print_combined_from_names(filename,
                                  core_xyz_filename,
                                  atoms_to_del,
                                  [bytes(name, encoding='utf-8') for name in frag_names])
    else:
        print_combined_from_xyz_filenames(filename,
                                          core_xyz_filename,
                                          atoms_to_del,
                                          [bytes(fn, encoding='utf-8') for fn in frag_xyz_filenames])
    return None
