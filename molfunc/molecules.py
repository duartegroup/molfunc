import numpy as np
from scipy.spatial import distance_matrix
from scipy.optimize import minimize
import networkx as nx
from molfunc import energy
from molfunc.atoms import smiles_to_atoms
from molfunc.atoms import xyz_file_to_atoms
from molfunc.bonds import get_avg_bond_length
from molfunc.atoms import NNAtom, Atom
from molfunc.geom import get_rotated_coords
from molfunc.exceptions import *
from molfunc.utils import requires_atoms
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')


class Molecule:

    def make_graph(self, rel_tolerance=0.2):
        """
        Make the molecular graph from the 'bonds' determined on a distance
        criteria. No distinction is made between single, double etc. bond types

        :param rel_tolerance: (float) Relative tolerance on what is classed
                                      as a bond. If a distance is
                                      < (1 + rel_tolerance) * r_avg
                                     then they are 'bonded'
        :return: None
        """

        graph = nx.Graph()
        for i in range(self.n_atoms):
            graph.add_node(i, atom_label=self.atoms[i].label)

        coordinates = self.get_coordinates()
        dist_mat = distance_matrix(coordinates, coordinates)

        # Loop over the unique pairs of atoms and add 'bonds'
        for i in range(self.n_atoms):
            for j in range(i+1, self.n_atoms):

                avg_bond_length = get_avg_bond_length(self.atoms[i].label,
                                                      self.atoms[j].label)

                # If the atoms are close enough add a bond (edge)
                if dist_mat[i, j] <= avg_bond_length * (1.0 + rel_tolerance):
                    graph.add_edge(i, j)

        self.graph = graph

        return None

    def set_atomic_valancies(self):
        """Set the atomic valency for each atom. Double/triple bonds are *not*
        distinct from single bonds"""

        for i in range(self.n_atoms):
            self.atoms[i].valence = len(list(self.graph.neighbors(i)))

        return None

    @requires_atoms()
    def translate(self, vec):
        """Translate the molecule by vector (np.ndarray, length 3)"""
        assert vec.shape == (3,)

        for atom in self.atoms:
            atom.translate(vec)
        return None

    @requires_atoms()
    def print_xyz_file(self):
        """Print a standard .xyz file from the Molecule's atoms"""

        if self.atoms is None or len(self.atoms) == 0:
            raise NoAtomsToPrint

        with open(f'{self.name}.xyz', 'w') as xyz_file:
            print(self.n_atoms, '\n', file=xyz_file)
            for atom in self.atoms:
                x, y, z = atom.coord
                print(f'{atom.label:<3}{x:^10.5f}{y:^10.5f}{z:^10.5f}',
                      file=xyz_file)

        return None

    @requires_atoms()
    def get_coordinates(self):
        """Return a n_atoms x 3 shape np.ndarray containing xyz coordinates"""
        return np.array([atom.coord for atom in self.atoms])

    def set_atoms(self, atoms):
        """Set the atoms (list(molfunc.atoms.Atom)) and the number of atoms"""
        assert type(atoms) is list
        if len(atoms) > 0:
            assert isinstance(atoms[0], Atom)

        self.atoms = atoms
        self.n_atoms = len(atoms)

        return None

    def __init__(self, name='mol', xyz_filename=None, smiles=None, atoms=None):
        """
        Base molecule class. Initialised in order of priority: SMILES string,
        xyz file, atoms

        ------------------------ Keyword Arguments ----------------------------
        :param name: (str)

        :param xyz_filename: (str) .xyz filename (or filepath) from which atoms
                             will be extracted

        :param smiles: (str) SMILES string defining the molecule from which a
                       3D structure as atoms are extracted using RDKit

        :param atoms: (list(molfunc.atom.Atom)) List of atoms used to
                      initialise the molecule
        """
        self.name = str(name)
        self.n_atoms = 0
        self.graph = None
        self.atoms = None

        if smiles is not None:
            # Use RDKit to convert SMILES -> atoms
            self.set_atoms(atoms=smiles_to_atoms(smiles))

        if xyz_filename is not None:
            # Initialisation with an xyz file takes precedence over SMILES
            self.set_atoms(atoms=xyz_file_to_atoms(xyz_filename))

        if atoms is not None:
            self.set_atoms(atoms)

        if self.n_atoms != 0:
            # If there are atoms in the molecule set the graph and valancies
            self.make_graph()
            self.set_atomic_valancies()


class CoreMolecule(Molecule):

    def _check_datom_idxs(self):
        """Ensure that all atoms to be replaced by fragments are monovalent"""

        for i in self.datom_idxs:

            if not 0 <= i < self.n_atoms:
                raise DatomsNotValid(f'Can\'t functionalise an atom {i} - not '
                                     f'in the list of atoms')

            if self.atoms[i].valence == 1:
                continue

            raise DatomsNotValid(f'Cannot modify atom {self.atoms[i].label} '
                                 f'with valency {self.atoms[i].valence}')
        return None

    def get_datom_nn(self, datom_idx):
        """
        Return the nearest neighbour atom to a particular atom to delete
        (datom) along with the shift vector e.g. for a datom_idx = 1, the
        nearest neighbour is C and the vector

            vec -->

           j     i                        atoms:
           C -- H                               C, 0, 0, 0
          /                                     H, 1, 0, 0
        H                                       H, -1, 0, -1


        :param datom_idx: (int) index of the atom to delete
        :return: (molfunc.atoms.NNatom)
        """

        for (i, j) in self.graph.edges:

            vec = self.atoms[i].coord - self.atoms[j].coord

            if i == datom_idx:
                return NNAtom(atom=self.atoms[j], shift_vec=vec)
            if j == datom_idx:
                return NNAtom(atom=self.atoms[i], shift_vec=-vec)

        return None

    def _delete_datoms(self):
        """Remove all datoms from the atoms list and set the atoms"""
        return self.set_atoms(atoms=[atom for i, atom in enumerate(self.atoms)
                                     if i not in self.datom_idxs])

    def __init__(self, name='mol', xyz_filename=None, smiles=None,
                 atoms_to_del=None, atoms=None):
        """
        Core molecule class

        :param name: (str)
        :param atoms_to_del: (list(int)) List of atom indexes to delete and
                             swap for a fragment. *Indexed from 1*
        :param xyz_filename: (str)
        :param smiles: (str) SMILES string
        :param atoms: (list(molfunc.atom.Atom)) List of atoms
        """
        super().__init__(name, xyz_filename, smiles, atoms)

        # Atom indexes to delete are indexed from 1. molfunc however indexes
        # atoms from 0 so atoms_to_del are the indexes minus 1
        self.datom_idxs = [i - 1 for i in atoms_to_del] if atoms_to_del is not None else []
        self._check_datom_idxs()

        # Nearest neighbour atoms to those deleted to enable translation of the
        # fragment
        self.nn_atoms = [self.get_datom_nn(i) for i in self.datom_idxs]

        # Remove the atoms in the datom_idxs list from the atoms
        self._delete_datoms()


class FragmentMolecule(Molecule):

    def get_ratom_nearest_neighbour(self):
        """
        Return the nearest neighbour atom to the atom with label 'R' the place
        holder atom that is deleted when the Fragment and Core molecules are
        bonded. It needs to be monovalent..

        :return: (molfunc.atoms.NNatom)
        """

        for (atom_i, atom_j) in self.graph.edges:

            if self.atoms[atom_i].label == 'R':
                if self.atoms[atom_i].valence != 1:
                    raise RAtomInvalidValence

                return NNAtom(atom=self.atoms[atom_j])

            if self.atoms[atom_j].label == 'R':
                if self.atoms[atom_j].valence != 1:
                    raise RAtomInvalidValence

                return NNAtom(atom=self.atoms[atom_i])

        raise RAtomNotFound

    def minimise_repulsion(self, other_mol, n=100, tolerance=0.1):
        """
        Minimise the 'energy' with respect to rigid body rotation of this
        molecule given some other (core) molecule

        :param other_mol: (molfunc.CoreMolecule)
        :param n: (int) Number of minimisation to perform as to (hopefully)
                  locate the global minimum
        :param tolerance: (float) Threshold on the energy minimisation 0.1 is
                          a reasonable compromise between speed and accuracy
        :return:
        """
        # Coords where the nearest neighbour atom is centered at the origin
        other_coords = other_mol.get_coordinates() - self.nn_atom.coord
        coords = self.get_coordinates() - self.nn_atom.coord

        # Set up initial parameters
        min_energy = best_x = None

        for _ in range(n):

            x0 = np.random.uniform(0.0, 2*np.pi, size=4)

            # Minimise repulsion between coords and the other_coords
            # TODO analytic derivative
            res = minimize(energy.repulsion, x0=x0, args=(coords, other_coords),
                           method='L-BFGS-B', tol=tolerance)

            if min_energy is None or res.fun < min_energy:
                min_energy = res.fun
                best_x = res.x

        # Set the new coordinates, with the best x array (lowest repulsion)
        new_coords = get_rotated_coords(x=best_x, coords=coords)
        for i in range(self.n_atoms):
            self.atoms[i].coord = new_coords[i] + self.nn_atom.coord

        return None

    def _delete_r_atom(self):
        """Delete the atom with label 'R' from the atoms"""
        return self.set_atoms(atoms=[atom for atom in self.atoms if atom.label != 'R'])

    def shift_to_core_atom(self, core_atom):
        """Given a core atom (molfunc.atoms.NNAtom) shift the fragment so the
        R atom nearest neighbour is the ideal length away from the core atom"""

        # Translate so the fragment nn atom is on top of the core atom
        self.translate(vec=core_atom.coord - self.nn_atom.coord)

        # Translate so the fragment has is the correct bond distance away
        ideal_bond_length = get_avg_bond_length(atom_i_label=core_atom.label,
                                                atom_j_label=self.nn_atom.label)
        self.translate(vec=core_atom.shift_vec * ideal_bond_length)

        # Update nn_atom coord as it will be used as the origin for rotation
        self.nn_atom.coord = core_atom.coord + core_atom.shift_vec * ideal_bond_length

        return None

    def __init__(self, name='mol', xyz_filename=None, smiles=None,
                 core_atom=None, atoms=None):
        """
        Fragment molecule class

        e.g. for a methyl self.atoms should be:
            C, 0, 0, 0
            R, 1, 0, 1
            H -1, 0, 1
            H 0, 1, -1
            H 0, -1, -1

        where the R atom will be removed and replaced for the core molecule.
        The C atom will be the nn_atom (closest to R)

        :param name: (str)
        :param core_atom: (molfunc.atoms.NNAtom) The atom on the core molecule
                          to which this fragment will be attached
        :param xyz_filename: (str)
        :param smiles: (str) SMILES string
        :param atoms: (list(molfunc.atom.Atom)) List of atoms
        """
        super().__init__(name, xyz_filename, smiles, atoms)

        self.nn_atom = self.get_ratom_nearest_neighbour()
        self._delete_r_atom()

        # To manipulate as a fragment there needs to be a core atom to
        # translate to
        if core_atom is not None:
            self.shift_to_core_atom(core_atom)


class CombinedMolecule(Molecule):

    def build(self):
        """Build the combined molecule by iterating through self.fragments
        fragments minimising the repulsion to the core mol at each point"""

        if len(self.fragments) != self.n_fragments_to_add:
            raise CombinationFailed('Number of fragments is not equal to the'
                                    'number of fragments to add (core datoms)')

        atoms = self.core_mol.atoms.copy()

        for i, fragment_mol in enumerate(self.fragments):
            fragment_mol.shift_to_core_atom(self.core_mol.nn_atoms[i])
            fragment_mol.minimise_repulsion(other_mol=self.core_mol)

            atoms += fragment_mol.atoms

        self.set_atoms(atoms)
        return None

    def __init__(self, core_mol, frag_smiles=None, frag_smiles_list=None,
                 name='mol', fragment=None, fragments=None):
        """
        Combined molecule class

        Fragments can be added from SMILES strings (e.g. C[*]), a list of
        SMILES strings (if the core atom is functionalised more than once),
        a FragmentMolecule or a list of FragmentMolecules again if there is
        more than 1 functionalisation to be performed

        :param name: (str)
        :param core_mol: (molfunc.molecules.CoreMolecule)
        :param frag_smiles: (str)
        :param frag_smiles_list: (list(str))
        :param fragment: (molfunc.molecules.FragmentMolecule)
        :param fragments: (list(molfunc.molecules.FragmentMolecule))
        """
        super().__init__(name=name)

        self.core_mol = core_mol
        self.n_fragments_to_add = len(core_mol.datom_idxs)

        if self.n_fragments_to_add == 0:
            raise CombinationFailed('Core molecule had no datoms')

        self.fragments = []
        if frag_smiles is not None:
            self.fragments = self.n_fragments_to_add * [FragmentMolecule(smiles=frag_smiles)]

        if frag_smiles_list is not None:
            self.fragments = [FragmentMolecule(smiles=smiles) for smiles in frag_smiles_list]

        if fragment is not None:
            assert isinstance(fragment, FragmentMolecule)
            self.fragments = self.n_fragments_to_add * [fragment]

        if fragments is not None:
            assert all(isinstance(fr, FragmentMolecule) for fr in fragments)
            self.fragments = fragments

        # If there are some fragments then build the combined molecule
        if len(self.fragments) > 0:
            self.build()
