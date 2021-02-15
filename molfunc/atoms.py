import numpy as np
import os
from molfunc.exceptions import *


class Atom:

    def __repr__(self):
        x, y, z = self.coord
        return f'[{self.label}, {x:.3f}, {y:.3f}, {z:.3f}]'

    def translate(self, vec):
        """Shift the atom in 3D space by a vector (np.ndarray) length 3"""
        self.coord += vec

    def __init__(self, atomic_symbol, x=0.0, y=0.0, z=0.0, coord=None):
        """
        Atom positioned in 3D space

        :param atomic_symbol: (str) Atomic symbol of this atom e.g. 'C' or 'Pd'
        :param x: (float) x coordinate of the atom
        :param y: (float)
        :param z: (float)

        :param coord: (np.ndarray)
        """
        self.label = str(atomic_symbol)

        if coord is not None:
            assert coord.shape == (3,)
            self.coord = coord

        else:
            self.coord = np.array([float(x), float(y), float(z)])

        # Atomic valence e.g. 4 for carbon in methane
        self.valence = 0


class NNAtom(Atom):

    def __init__(self, atom, shift_vec=None):
        """
        Nearest neighbour atom

        :param atom: (molfunc.atoms.Atom)
        :param shift_vec: (np.ndarray) length 3
        """
        super().__init__(atomic_symbol=atom.label, coord=atom.coord)
        self.shift_vec = shift_vec

        if self.shift_vec is not None:
            # Ensure the vector is normalised
            self.shift_vec = shift_vec / np.linalg.norm(shift_vec)


def xyz_file_to_atoms(filename):
    """
    From an .xyz file get a list of atoms. If a fragment set of atoms they
    should either contain an R atom which will be deleted and used to generate
    the functionalised molecule

    :param filename: (str)
    :return: (list(molfunc.atoms.Atom))
    """
    atoms = []

    if not os.path.exists(filename):
        raise XYZfileDidNotExist

    if not filename.endswith('.xyz'):
        raise NotXYZfile

    with open(filename, 'r') as xyz_file:
        try:
            # First item in an xyz file is the number of atoms
            n_atoms = int(xyz_file.readline().split()[0])

            # XYZ lines should be the following 2 + n_atoms lines
            xyz_lines = xyz_file.readlines()[1:n_atoms+1]

        except (IndexError, ValueError):
            raise XYZfileMalformatted

        for line in xyz_lines:

            try:
                atom_label, x, y, z = line.split()[:4]
                atoms.append(Atom(atomic_symbol=atom_label, x=x, y=y, z=z))

            except (IndexError, ValueError):
                raise XYZfileMalformatted

    return atoms


def smiles_to_atoms(smiles, n_confs=1):
    """
    From a SMILES string generate n_confs conformers using RDKit with the
    ETKDGv2 algorithm

    If there is an * atom in the SMILES string this function will convert it
    to a Li atom to generate the structure with RDKit then swap it for R from
    which will be added to a core

    :param smiles: (str)
    :param n_confs: (int)
    :return: (list(list(molfunc.atoms.Atom)) if n_confs > 1 else
             (list(molfunc.atoms.Atom))
    """
    # Allow for molfunc usage without RDKit if no SMILES -> 3D conversion is
    # required
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Swap Fr for Li so RDKit generates something reasonable
    if 'Li' in smiles:
        raise MolFuncCritical('Cannot convert a molecule containing Li')

    if '[*]' in smiles:
        smiles = smiles.replace('[*]', '[Li]')

    # Generate an RDKit Molecule object, add hydrogens and generate conformers
    try:
        mol_obj = Chem.MolFromSmiles(smiles)
        mol_obj = Chem.AddHs(mol_obj)
        method = AllChem.ETKDGv2()
        method.pruneRmsThresh = 0.3
        conf_ids = list(AllChem.EmbedMultipleConfs(mol_obj, numConfs=n_confs,
                                                   params=method))
    except:
        # RDKit can raise unhelpful exceptions
        raise RDKitFailed

    # If n_confs > 1 then append to this list the list of Atoms
    atoms_list = []

    for i in range(len(conf_ids)):
        mol_block_lines = Chem.MolToMolBlock(mol_obj,
                                             confId=conf_ids[i]).split('\n')
        atoms = []

        # Get the atomic symbols and atomic coordinates from the block
        for line in mol_block_lines:
            items = line.split()
            if len(items) == 16:
                atom_label, x, y, z = items[3], items[0], items[1], items[2]

                # Swap Li to R
                if atom_label == 'Li':
                    atom_label = 'R'

                atoms.append(Atom(atomic_symbol=atom_label, x=x, y=y, z=z))

        # If only one conformer is requested then only return a list of atoms
        if n_confs == 1:
            return atoms

        atoms_list.append(atoms)

    return atoms_list


# vdw radii from https://books.google.no/books?id=bNDMBQAAQBAJ
vdw_radii = {'H': 1.1, 'He': 1.4, 'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.7,
             'N': 1.55, 'O': 1.52, 'F': 1.47, 'Ne': 1.54, 'Na': 2.27,
             'Mg': 1.73, 'Al': 1.84, 'Si': 2.1, 'P': 1.8, 'S': 1.8,
             'Cl': 1.75, 'Ar': 1.88, 'K': 2.75, 'Ca': 2.31, 'Sc': 2.15,
             'Ti': 2.11, 'V': 2.07, 'Cr': 2.06, 'Mn': 2.05, 'Fe': 2.04,
             'Co': 2.0, 'Ni': 1.97, 'Cu': 1.96, 'Zn': 2.01, 'Ga': 1.87,
             'Ge': 2.11, 'As': 1.85, 'Se': 1.9, 'Br': 1.85, 'Kr': 2.02,
             'Rb': 3.03, 'Sr': 2.49, 'Y': 2.32, 'Zr': 2.23, 'Nb': 2.18,
             'Mo': 2.17, 'Tc': 2.16, 'Ru': 2.13, 'Rh': 2.1, 'Pd': 2.1,
             'Ag': 2.11, 'Cd': 2.18, 'In': 1.93, 'Sn': 2.17, 'Sb': 2.06,
             'Te': 2.06, 'I': 1.98, 'Xe': 2.16, 'Cs': 3.43, 'Ba': 2.68,
             'La': 2.43, 'Ce': 2.42, 'Pr': 2.4, 'Nd': 2.39, 'Pm': 2.38,
             'Sm': 2.36, 'Eu': 2.35, 'Gd': 2.34, 'Tb': 2.33, 'Dy': 2.31,
             'Ho': 2.3, 'Er': 2.29, 'Tm': 2.27, 'Yb': 2.26, 'Lu': 2.24,
             'Hf': 2.23, 'Ta': 2.22, 'W': 2.18, 'Re': 2.16, 'Os': 2.16,
             'Ir': 2.13, 'Pt': 2.13, 'Au': 2.14, 'Hg': 2.23, 'Tl': 1.96,
             'Pb': 2.02, 'Bi': 2.07, 'Po': 1.97, 'At': 2.02, 'Rn': 2.2,
             'Fr': 3.48, 'Ra': 2.83, 'Ac': 2.47, 'Th': 2.45, 'Pa': 2.43,
             'U': 2.41, 'Np': 2.39, 'Pu': 2.43, 'Am': 2.44, 'Cm': 2.45,
             'Bk': 2.44, 'Cf': 2.45, 'Es': 2.45, 'Fm': 2.45, 'Md': 2.46,
             'No': 2.46, 'Lr': 2.46}
