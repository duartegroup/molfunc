#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <numeric>
#include "memory"
#include "molecules.h"
#include "utils.h"


using namespace std;

namespace molfunc{

    Molecule::Molecule() = default;

    Molecule::Molecule(const string& xyz_filename) {
        /*********************************************************
         * Construct a Molecule from a standard .xyz file
         * https://en.wikipedia.org/wiki/XYZ_file_format.
         * No requirement for charge or spin state attributes,
         * as molfunc only operates on geometries
         *
         * Arguments:
         *      xyz_filename (string):
         *
         * Example:
         *      Molecule mol = Molecule("methane.xyz");
         *      // mol.n_atoms() -> 5
         ********************************************************/

        if (!utils::ends_with(xyz_filename, ".xyz")){
            throw runtime_error("Expecting a .xyz file, had: "+xyz_filename);
        }

        set_atoms(xyz_filename);
        assign_coordinates();
        construct_graph();
    }

    Molecule::Molecule(const vector<Atom>& atoms){
        /*********************************************************
         * Construct a Molecule from a set of atoms
         *
         * Arguments:
         *      Atom (list(Atom)):
         *
         * Example:
         *      vector<Atom> atoms = {Atom("H", 0.0, 0.0, 0.0)};
         *      Molecule mol = Molecule(atoms);
         *      // mol.n_atoms() -> 1
         ********************************************************/

        this->atoms = atoms;
        assign_coordinates();
        construct_graph();
    }

    void Molecule::set_atoms(const string& xyz_filename){
        /*********************************************************
         * Construct a Molecule from a standard .xyz file
         * https://en.wikipedia.org/wiki/XYZ_file_format
         *
         * Arguments:
         *      xyz_filename (string):
         ********************************************************/

        string line, item;
        int line_n = 0;
        ifstream xyz_file (xyz_filename);

        if (!xyz_file.is_open()){
            string file_path = __FILE__;
            string dir_path = file_path.substr(0, file_path.rfind('/'));
            throw runtime_error("Failed to open: "+xyz_filename);
        }

        int decl_n_atoms = 0;  // Declared number of atoms

        // Iterate through the xyz file
        while (getline(xyz_file, line, '\n')) {
            line_n += 1;


            if (line_n == 1) {          // Number of atoms line
                decl_n_atoms = stoi(line);
                continue;
            }

            if (line_n == 2){           // Title line
                xyz_title_line = line;
                continue;
            }

            if (line.empty()) continue; // Skip any blank lines

            vector<string> xyz_items = utils::split(line, ' ');
            if (xyz_items.size() != 4){
                throw runtime_error("Malformatted xyz file line: " + line);
            }

            atoms.emplace_back(xyz_items[0],         // Atomic symbol
                               stod(xyz_items[1]),   // x
                               stod(xyz_items[2]),   // y
                               stod(xyz_items[3]));  // z

        }
        xyz_file.close();

        if (decl_n_atoms != n_atoms()) {
            throw runtime_error("Number of atoms "+to_string(n_atoms())+
            " not equal to the number declared "
            + to_string(decl_n_atoms));
        }


    }

    void Molecule::assign_coordinates(){
        /*******************************************************************
         * Assign the coordinates for this molecule in a memory-contiguous
         * way
         ******************************************************************/
        coordinates.reserve(n_atoms());

        for (auto &atom : atoms){
            coordinates.push_back({atom.x(), atom.y(), atom.z()});

            atom.ptr_x = &coordinates.back()[0];
            atom.ptr_y = &coordinates.back()[1];
            atom.ptr_z = &coordinates.back()[2];
        }
    }

    double Molecule::distance(unsigned long i, unsigned long j){
        /*******************************************************************
         * Distance between two atoms
         *
         * Arguments:
         *      i (int): Atom index
         *
         *      j (int): Atom index
         *
         * Example:
         *      Molecule mol = Molecule("methane.xyz");
         *      // Methane molecule with the carbon as the first atom index
         *      // mol.distance(0, 1) -> ~1.0
         ******************************************************************/

        if (i >= n_atoms() ||j >= n_atoms() || coordinates.empty()){
            throw runtime_error("Invalid index: "+to_string(i)+" or "+ to_string(j)+
                                " must be present in the molecule");
        }

        double sq_dist = 0.0;

        for (int k=0; k<3; k++){
            double tmp = coordinates[i][k] - coordinates[j][k];
            sq_dist += tmp * tmp;
        }

        return sqrt(sq_dist);
    }

    bool Molecule::is_bonded_on_distance(unsigned long i, unsigned long j) {
        /*******************************************************************
         * Are two atoms bonded based on a distance criteria?
         * Uses covalent radii from: Dalton Trans., 2008, 2832
         *
         *
         * Arguments:
         *      i (int): Atom index
         *
         *      j (int): Atom index
         *
         * Example:
         *      Molecule mol = Molecule("methane.xyz");
         *      // Methane molecule with the carbon as the first atom index
         *      // mol.is_bonded_on_distance(0, 1) -> true
         ******************************************************************/

        //                           H      He    B   ...
        vector<double> cov_radii = {0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11,
                                    1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.61, 1.52, 1.50, 1.24,
                                    1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54,
                                    1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15,
                                    2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.95, 1.94, 1.92, 1.92, 1.89, 1.90,
                                    1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48,
                                    1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69};
        //                                                                                                  ^ Cm

        if (i == j){
            throw runtime_error("Self bonded"+to_string(i)+"-"+to_string(i)+
                                " have no bond definition");
        }

        if (i >= n_atoms() || j >= n_atoms()){
            throw runtime_error("Cannot determine if "+
                                to_string(i)+"-"+ to_string(j)+" are bonded."
                                " At least one was outside the set of atoms");
        }

        if (atoms[i].is_dummy() || atoms[j].is_dummy()){
            // Dummy atoms are by definition bonded to the closest atom
            return false;
        }

        unsigned int z_i = atoms[i].atomic_number();
        unsigned int z_j = atoms[j].atomic_number();

        if (z_i > cov_radii.size() || z_j > cov_radii.size()){
            throw runtime_error("Unknown element "+atoms[i].symbol+" "
                                "or "+atoms[j].symbol+". Not in known covalent radii");
        }

        // Is the distance less than 1.2x the sum of the covalent radii?
        // (where the index is one minus the atomic number (indexed from 1))
        return distance(i, j) < 1.2*(cov_radii[z_i-1] + cov_radii[z_j-1]);
    }

    void Molecule::construct_graph() {
        /*********************************************************
         * Using a set of atoms defined for this molecule set the
         * nodes as atoms (with the correct indexing) and edges
         * as 'covalent' bonds between them
         ********************************************************/

        graph = Graph();
        for (unsigned long i=0; i < n_atoms(); i++){
            graph.add_node(Node(i, atoms[i].symbol));
        }

        // Enumerate all unique pairwise atoms and determine if they
        // are 'bonded'
        for (unsigned long i=0; i < n_atoms(); i++){
            for (unsigned long j=i+1; j < n_atoms(); j++){
                if (is_bonded_on_distance(i, j)){
                    graph.add_edge(i, j);
                }
            }// j
        } // i

        // Bond any dummy atoms to the closest neighbour
        double min_distance = INFINITY;
        unsigned long closest_idx;

        for (unsigned long i=0; i < n_atoms(); i++){

            if (!atoms[i].is_dummy()) continue;

            for (unsigned long j=0; j < n_atoms(); j++){

                double dist = distance(i, j);

                if (i != j && dist < min_distance) {
                    closest_idx = j;
                    min_distance = dist;
                }
            }// j

            graph.add_edge(i, closest_idx);
        } // i
    }

    unsigned long Molecule::n_atoms() const {
        // Number of atoms in this molecule
        return atoms.size();
    }

    unsigned long Molecule::n_masked_atoms(){
        // Number of masked atoms in this system

        unsigned long n_atoms = 0;
        for (auto &atom : atoms){
            if (atom.masked) n_atoms++;
        }
        return n_atoms;
    }

    unsigned long Molecule::n_unmasked_atoms() {return n_atoms() - n_masked_atoms();}

    void Molecule::print_xyz_file(const string& filename){
        /*********************************************************
         * Generate a standard .xyz file for a molecule
         *
         * Arguments:
         *      filename (str):
         ********************************************************/

        if (atoms.empty()){
            throw runtime_error("Could not print a .xyz file- had no atoms");
        }

        ofstream xyz_file (filename);

        if (xyz_file.is_open()){

            xyz_file << fixed;
            xyz_file << setprecision(6);

            // ---------------------------------------------------
            xyz_file << to_string(n_unmasked_atoms()) << '\n'
                     << "molfunc generated" << '\n';

            for (auto &atom: atoms){

                if (atom.masked) continue;  // Skip masked atoms

                xyz_file << atom.symbol   << "    "
                         << atom.x() << "    "
                         << atom.y() << "    "
                         << atom.z() << "    " <<  '\n';
            }
            // ---------------------------------------------------

            xyz_file.close();
        }

        else throw runtime_error("Cannot open "+filename);
    }

    CoreMolecule::CoreMolecule(const string &xyz_filename):
                              Molecule(xyz_filename){
        /*********************************************************
         * Construct a CoreMolecule from a standard .xyz file that
         * contains one or more "R" dummy atoms that will be
         * replaced by a fragment(s)
         *
         * Arguments:
         *      xyz_filename (string):
         ********************************************************/
         if (n_atoms() == n_unmasked_atoms()){
             throw runtime_error("Cannot construct a CoreMolecule with no "
                                 "dummy (R) atoms present in the xyz file");
         }

    }

    CoreMolecule::CoreMolecule(const string &xyz_filename,
                               const vector<unsigned int>& atoms_to_del)
                               : Molecule(xyz_filename) {
        /*********************************************************
         * Construct a CoreMolecule from a standard .xyz file and a
         * number of atoms to delete, which will be replaced by
         * fragments.
         *
         * Arguments:
         *      xyz_filename (string):
         *
         *      atoms_to_del (list(int)): Atom indexes
         *                                (indexed from 0)
         ********************************************************/

        for (auto atom_idx : atoms_to_del){

            if (atom_idx >= n_atoms()){
                throw out_of_range("Cannot delete atom " + to_string(atom_idx) +
                                   " . Not present in " + xyz_filename);
            }

            if (graph.n_neighbours(atom_idx) != 1){
                throw runtime_error("Deleted atoms must be monovalent. Atom "+
                                    to_string(atom_idx)+" was not.");
            }

            atoms[atom_idx].masked = true;
        }
    }
}

