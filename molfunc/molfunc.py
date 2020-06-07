import argparse
from molfunc.molecules import CoreMolecule
from molfunc.molecules import CombinedMolecule


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('xyz_filename', action='store', type=str,
                        help='<Required> .xyz file to functionalise. '
                             'e.g. methane.xyz')
    parser.add_argument('-a', '--atom_ids', action='store', nargs='+', type=str,
                        help='<Required> List of atom ids to swap for '
                             'fragments. e.g. 1 2')
    parser.add_argument('-s', '--fragment_smiles', action='store', type=str,
                        help='<Required> SMILES string of the fragment to add.'
                             ' e.g. [*]C')

    return parser.parse_args()


def main():

    args = get_args()
    atom_ids = [int(i) for i in args.atom_ids]

    core_mol = CoreMolecule(xyz_filename=args.xyz_filename,
                            atoms_to_del=atom_ids)

    joined = CombinedMolecule(core_mol=core_mol,
                              frag_smiles=args.fragment_smiles,
                              name=args.xyz_filename.replace('.xyz', '_mod'))
    joined.print_xyz_file()
