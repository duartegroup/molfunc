import numpy as np
from scipy.spatial import distance_matrix
from scipy.optimize import minimize
import networkx as nx
from ffunc.atoms import smiles_to_atoms
from ffunc.atoms import xyz_file_to_atoms
from ffunc.bonds import get_avg_bond_length
from ffunc.atoms import NNAtom
from ffunc.geom import rotation_matrix
from ffunc.exceptions import *
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')


def get_rotated_coords(x, theta, coords):
    rot_mat = rotation_matrix(axis=x/np.linalg.norm(x), theta=theta)
    return [np.matmul(rot_mat, coord) for coord in coords]


def energy_func(x, theta, coords, alt_coords):
    """Evaluate the energy for coords that have been rotated by theta radians in the axis x"""

    coords = get_rotated_coords(x, theta, coords)
    dist_mat = distance_matrix(coords, alt_coords)
    energy = np.sum(np.power(dist_mat, -4))

    return energy


class Molecule:

    def _make_graph(self, rel_tolerance=0.2):
        """
        Make the molecular graph from the 'bonds' determined on a distance criteria. No distinction is made between
        single, double etc. bond types

        :param rel_tolerance: (float)
        :return: None
        """

        graph = nx.Graph()
        for i in range(self.n_atoms):
            graph.add_node(i, atom_label=self.atoms[i].label)

        coordinates = self.get_coords()
        dist_mat = distance_matrix(coordinates, coordinates)

        # Loop over the unique pairs of atoms and add 'bonds'
        for i in range(self.n_atoms):
            for j in range(i+1, self.n_atoms):

                avg_bond_length = get_avg_bond_length(atom_i_label=self.atoms[i].label,
                                                      atom_j_label=self.atoms[j].label)

                # If the distance between atoms i and j are less or equal to 1.2x average length add a 'bond'
                if dist_mat[i, j] <= avg_bond_length * (1.0 + rel_tolerance):
                    graph.add_edge(i, j)

        self.graph = graph

        return None

    def _set_atomic_valancies(self):
        """Set the atomic valency for each atom. Double/triple bonds are *not* distinct from single bonds"""

        for (atom_i, atom_j) in self.graph.edges:

            self.atoms[atom_i].valence += 1
            self.atoms[atom_j].valence += 1

        return None

    def translate(self, vec):
        """Translate the molecule by vector (np.ndarray, length 3)"""
        for atom in self.atoms:
            atom.translate(vec)
        return None

    def print_xyz_file(self):
        """Print a standard xyz file from the Molecule's atoms"""

        if self.atoms is None or len(self.atoms) == 0:
            raise NoAtomsToPrint

        with open(f'{self.name}.xyz', 'w') as xyz_file:
            print(self.n_atoms, '\n', file=xyz_file)
            for atom in self.atoms:
                x ,y, z = atom.coord
                print(f'{atom.label:<3}{x:^10.5f}{y:^10.5f}{z:^10.5f}', file=xyz_file)

        return None

    def get_coords(self):
        """Return a np.ndarray of size n_atoms x 3 containing the xyz coordinates of the molecule"""
        return np.array([atom.coord for atom in self.atoms])

    def set_atoms(self, atoms):

        self.atoms = atoms
        self.n_atoms = len(atoms)

        return None

    def __init__(self, name='molecule', xyz_filename=None, smiles=None, atoms=None):

        self.name = name
        self.n_atoms = 0
        self.graph = None
        self.atoms = None

        if smiles is not None:
            self.set_atoms(atoms=smiles_to_atoms(smiles))

        if xyz_filename is not None:
            # Initialisation with an xyz file takes precedence over SMILES string
            self.set_atoms(atoms=xyz_file_to_atoms(xyz_filename))

        if atoms is not None:
            self.set_atoms(atoms)

        if self.n_atoms == 0:
            # Cannot continue if there are no atoms in the molecule..
            return

        self._make_graph()
        self._set_atomic_valancies()


class CoreMolecule(Molecule):

    def _check_datom_idxs(self):
        """Ensure that all the atoms that will be replaced by fragments are monovalent"""

        for i in self.datom_idxs:
            if self.atoms[i].valence > 1:
                exit(f'Cannot modify atom {self.atoms[i].label} with valency {self.atoms[i].valence }')

        return None

    def get_datom_nearest_neighbour(self, datom_idx):
        """
        Return the nearest neighbour atom to a particular atom to delete (datom)

        :param datom_idx: (int) index of the atom to delete
        :return: (ffunc.NNatom)
        """

        for (atom_i, atom_j) in self.graph.edges:

            vec = self.atoms[atom_i].coord - self.atoms[atom_j].coord

            if atom_i == datom_idx:
                return NNAtom(atom=self.atoms[atom_j], shift_vec=vec)
            if atom_j == datom_idx:
                return NNAtom(atom=self.atoms[atom_i], shift_vec=-vec)

        return None

    def _delete_atoms(self):
        """Delete all datoms from the atoms list"""
        return self.set_atoms(atoms=[atom for i, atom in enumerate(self.atoms) if i not in self.datom_idxs])

    def __init__(self, name='molecule', xyz_filename=None, smiles=None, atoms_to_del=None):

        super(CoreMolecule, self).__init__(name=name, xyz_filename=xyz_filename, smiles=smiles)

        # Atom indexes to delete are the atoms minus one as atoms_to_del should not have 0 in the list
        self.datom_idxs = [i - 1 for i in atoms_to_del] if atoms_to_del is not None else None
        self._check_datom_idxs()

        # Nearest neighbour atoms to those deleted to enable translation of the fragment
        self.nn_atoms = [self.get_datom_nearest_neighbour(datom_idx) for datom_idx in self.datom_idxs]
        self._delete_atoms()


class FragmentMolecule(Molecule):

    def get_ratom_nearest_neighbour(self):
        """
        Return the nearest neighbour atom to the atom with label 'R'

        :return: (ffunc.NNatom)
        """

        for (atom_i, atom_j) in self.graph.edges:

            if self.atoms[atom_i].label == 'R':
                if self.atoms[atom_i].valence > 1:
                    exit('R atoms must have valency 1')

                return NNAtom(atom=self.atoms[atom_j])

            if self.atoms[atom_j].label == 'R':
                if self.atoms[atom_j].valence > 1:
                    exit('R atoms must have valency 1')

                return NNAtom(atom=self.atoms[atom_i])

        raise RAtomNotFound

    def minimise_repulsion(self, other_mol, n=50):

        # Get the where the nearest neighbour atom is centered at the origin
        other_coords = other_mol.get_coords() - self.nn_atom.coord

        coords = self.get_coords() - self.nn_atom.coord

        min_energy = 9999999.9
        best_x, best_theta = None, None

        for _ in range(n):

            theta, x0 = np.random.uniform(0.0, 2*np.pi), np.random.uniform(0.0, 1.0, size=3)
            res = minimize(energy_func, x0=x0, args=(theta, coords, other_coords), method='L-BFGS-B')

            if res.fun < min_energy:
                # print('E = ', min_energy)
                min_energy = res.fun
                best_x, best_theta = res.x, theta

        new_coords = get_rotated_coords(x=best_x, theta=best_theta, coords=coords)
        for i in range(self.n_atoms):
            self.atoms[i].coord = new_coords[i] + self.nn_atom.coord

        return None

    def _delete_r_atom(self):
        """Delete the atom with label 'R' from the atoms"""
        return self.set_atoms(atoms=[atom for atom in self.atoms if atom.label != 'R'])

    def __init__(self, name='molecule', xyz_filename=None, smiles=None, core_atom=None):
        super(FragmentMolecule, self).__init__(name=name, xyz_filename=xyz_filename, smiles=smiles)

        # To manipulate as a fragment there needs to be a core atom to translate to
        if core_atom is None:
            return

        self.nn_atom = self.get_ratom_nearest_neighbour()
        self._delete_r_atom()

        # Translate so the fragment nn atom is on top of the core atom
        self.translate(vec=core_atom.coord - self.nn_atom.coord)

        # Translate so the fragment has is the correct bond distance away
        ideal_bond_length = get_avg_bond_length(atom_i_label=core_atom.label, atom_j_label=self.nn_atom.label)
        self.translate(vec=core_atom.shift_vec * ideal_bond_length)

        # Update the nn_atom coord as it will be used as the origin for rotation
        self.nn_atom.coord = core_atom.coord + core_atom.shift_vec * ideal_bond_length


class CombinedMolecule(Molecule):

    def _check_fragment_smiles_list(self):
        """The number of fragments smiles in the list *must* be equal to the number of deleted atoms"""
        if len(self.frag_smiles_list) != self.n_fragments_to_add:
            raise FragmentSMILESListMalformatted

    def __init__(self, core_mol, frag_smiles=None, frag_smiles_list=None, name='molecule'):

        super(CombinedMolecule, self).__init__(name=name)
        self.core_mol = core_mol
        self.n_fragments_to_add = len(self.core_mol.datom_idxs)

        if frag_smiles_list is not None:
            self.frag_smiles_list = frag_smiles_list
            self._check_fragment_smiles_list()

        elif frag_smiles is not None:
            # If the full list of fragment SMILES is not specified then populate with the number of atoms
            # that have been deleted from the CoreMolecule
            self.frag_smiles_list = self.n_fragments_to_add * [frag_smiles]

        else:
            exit('Could not create combined molecule')

        atoms = self.core_mol.atoms

        for i, fragment_smiles in enumerate(self.frag_smiles_list):
            fragment_mol = FragmentMolecule(smiles=fragment_smiles, core_atom=core_mol.nn_atoms[i])
            fragment_mol.minimise_repulsion(other_mol=self.core_mol)

            atoms += fragment_mol.atoms

        self.set_atoms(atoms)
