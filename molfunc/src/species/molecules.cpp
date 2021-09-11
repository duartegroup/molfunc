#include <fstream>
#include "string"
#include <cmath>
#include "memory"
#include "species/molecules.h"
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
        construct_graph();
    }

    Molecule::Molecule(vector<Atom3D>& atoms){
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

        this->atoms = {};

        for (auto &atom3d : atoms){
            this->coordinates.push_back({atom3d.x(), atom3d.y(), atom3d.z()});
            this->atoms.emplace_back(atom3d.symbol);
        }

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

        unsigned long decl_n_atoms = 0;  // Declared number of atoms

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

            atoms.emplace_back(xyz_items[0]);         // Atomic symbol
            coordinates.push_back({stod(xyz_items[1]),   // x
                                  stod(xyz_items[2]),    // y
                                  stod(xyz_items[3])});  // z
        }
        xyz_file.close();

        if (decl_n_atoms != n_atoms()) {
            throw runtime_error("Number of atoms "+to_string(n_atoms())+
            " not equal to the number declared "
            + to_string(decl_n_atoms));
        }


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

        // Is the distance less than 1.2x the sum of the covalent radii?
        auto max_dist = 1.2*(atoms[i].covalent_radius()
                             + atoms[j].covalent_radius());

        return distance(i, j) < max_dist;
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

    CoreMolecule::CoreMolecule() = default;

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

