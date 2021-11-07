import argparse
from molfunc.molecules import print_combined_molecule, print_all_combined_molecules


def get_args():
    """Get the command line arguments using argparse"""

    parser = argparse.ArgumentParser()
    parser.add_argument('xyz_filename',
                        action='store',
                        type=str,
                        help='<Required> .xyz file to functionalise. '
                             'e.g. methane.xyz')

    parser.add_argument('-a', '--atom_ids',
                        action='store',
                        nargs='+',
                        default=[],
                        type=str,
                        help='<Required> List of atom ids to swap for '
                             'fragments. e.g. 1 2')

    parser.add_argument('--all',
                        action='store_true',
                        default=False,
                        help='Attempt to build all combinations of '
                             'fragments for a defined core')

    parser.add_argument('-s', '--fragment_smiles',
                        action='store',
                        nargs='+',
                        type=str,
                        help='SMILES string of the fragment(s) to add.'
                             ' e.g. [*]C')

    parser.add_argument('-f', '--fragment_names',
                        action='store',
                        nargs='+',
                        type=str,
                        help='Names of the fragment(s) to add. e.g. Me')

    parser.add_argument('-n', '--name',
                        action='store',
                        type=str,
                        default=None,
                        help='Name of the combined molecule')

    return parser.parse_args()


def main():

    args = get_args()

    # Default name as the core name but modified ('mod')
    if args.name is None:
        args.name = args.xyz_filename.replace('.xyz', '_mod')

    if args.all:
        print_all_combined_molecules(core_xyz_filename=args.xyz_filename,
                                     atoms_to_del=[int(i) for i in args.atom_ids],
                                     name=args.name)

    else:
        print_combined_molecule(core_xyz_filename=args.xyz_filename,
                                atoms_to_del=[int(i) for i in args.atom_ids],
                                frag_smiles=args.fragment_smiles,
                                frag_names=args.fragment_names,
                                name=args.name)
