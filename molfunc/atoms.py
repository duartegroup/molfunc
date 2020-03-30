import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from molfunc.exceptions import *


class Atom:

    def __repr__(self):
        return f'[{self.label}, {self.coord[0]:.4f}, {self.coord[1]:.4f}, {self.coord[2]:.4f}]'

    def translate(self, vec):
        self.coord += vec

    def __init__(self, atomic_symbol, x, y, z):

        self.label = atomic_symbol
        self.coord = np.array([x, y, z])

        self.valence = 0


class NNAtom(Atom):

    def __init__(self, atom, shift_vec=None):
        super(NNAtom, self).__init__(atomic_symbol=atom.label, x=atom.coord[0], y=atom.coord[1],
                                     z=atom.coord[2])

        if shift_vec is not None:
            self.shift_vec = shift_vec / np.linalg.norm(shift_vec)      # Ensure the vector is normalised


def xyz_file_to_atoms(filename):
    """/
    From an .xyz file get a list of atoms

    :param filename: (str)
    :return: (list(molfunc.Atom))
    """

    atoms = []

    if not os.path.exists(filename):
        raise XYZfileDidNotExist

    if not filename.endswith('.xyz'):
        raise XYZfileWrongFormat

    # Open the file that exists and should(!) be in the correct format
    with open(filename, 'r') as xyz_file:
        n_atoms = int(xyz_file.readline().split()[0])       # First item in an xyz file is the number of atoms
        xyz_lines = xyz_file.readlines()[1:n_atoms+1]       # XYZ lines should be the following 2 + n_atoms lines

        for line in xyz_lines:

            try:
                atom_label, x, y, z = line.split()[:4]
                atoms.append(Atom(atomic_symbol=atom_label, x=float(x), y=float(y), z=float(z)))

            except (IndexError, TypeError):
                raise XYZfileMalformatted

    return atoms


def smiles_to_atoms(smiles, n_confs=1):
    """
    From a SMILES string generate n_confs conformers using RDKit

    :param smiles: (str)
    :param n_confs: (int)
    :return: (list(list(molfunc.Atom)) if n_confs > 1 else (list(fffunc.Atom))
    """

    # Generate an RDKit Molecule object, add hydrogens and generate conformers
    mol_obj = Chem.MolFromSmiles(smiles)
    mol_obj = Chem.AddHs(mol_obj)
    method = AllChem.ETKDGv2()
    method.pruneRmsThresh = 0.3
    conf_ids = list(AllChem.EmbedMultipleConfs(mol_obj, numConfs=n_confs, params=method))

    # If n_confs > 1 then generate a list of Atoms
    atoms_list = []

    for i in range(len(conf_ids)):
        mol_block_lines = Chem.MolToMolBlock(mol_obj, confId=conf_ids[i]).split('\n')
        atoms = []

        for line in mol_block_lines:
            split_line = line.split()
            if len(split_line) == 16:
                atom_label, x, y, z = split_line[3], split_line[0], split_line[1], split_line[2]

                atoms.append(Atom(atomic_symbol=atom_label, x=float(x), y=float(y), z=float(z)))

        # If only one conformer is requested then only return a list of atoms
        if n_confs == 1:
            return atoms

        atoms_list.append(atoms)

    return atoms_list
