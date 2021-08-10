#include <fstream>
#include <string>
#include <cmath>
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

        set_atoms(xyz_filename);
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
        ifstream xyz_file (xyz_filename);

        // Declared number of atoms
        int decl_n_atoms = 0;

        // Iterate through the xyz file
        while (getline(xyz_file, line, '\n')) {

            if (line.empty()) {   // Ignore any blank lines
                continue;
            }

            if (decl_n_atoms == 0) {
                decl_n_atoms = stoi(line);
                continue;
            }

            vector<string> xyz_items = utils::split(line, ' ');
            if (xyz_items.size() != 4){
                throw runtime_error("Malformatted xyz file, expecting a "
                                    "*A  0.0 0.0 0.0* structure to the line");
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

        if (i >= n_atoms() ||j >= n_atoms()){
            throw runtime_error("Invalid index: "+to_string(i)+" or "+ to_string(j)+
                                "must be present in the molecule");
        }

        double sq_dist = 0.0;

        for (int k=0; k<3; k++){
            double tmp = atoms[i].coord[k] - atoms[j].coord[k];
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
    }

    unsigned long Molecule::n_atoms() const {
        return atoms.size();
    }

}

