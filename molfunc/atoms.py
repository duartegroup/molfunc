import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
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
    From an .xyz file get a list of atoms

    :param filename: (str)
    :return: (list(molfunc.Atom))
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

    :param smiles: (str)
    :param n_confs: (int)
    :return: (list(list(molfunc.atoms.Atom)) if n_confs > 1 else
             (list(molfunc.atoms.Atom))
    """
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

                atoms.append(Atom(atomic_symbol=atom_label, x=x, y=y, z=z))

        # If only one conformer is requested then only return a list of atoms
        if n_confs == 1:
            return atoms

        atoms_list.append(atoms)

    return atoms_list
