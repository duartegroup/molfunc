#include "species.h"
#include <fstream>
#include <string>
#include <iomanip>
#include "cmath"
#include "memory"


namespace molfunc{

    Species::Species() = default;

    unsigned long Species::n_atoms() const {
        // Number of atoms in this molecule
        return atoms.size();
    }

    unsigned long Species::n_masked_atoms(){
        // Number of masked atoms in this system

        unsigned long n_atoms = 0;
        for (auto &atom : atoms){
            if (atom.masked) n_atoms++;
        }
        return n_atoms;
    }

    unsigned long Species::n_unmasked_atoms() {return n_atoms() - n_masked_atoms();}

    double Species::distance(unsigned long i, unsigned long j){
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

    void Species::translate(array<double, 3> vec){
        /*******************************************************************
         * Translate all the atoms in this molecule by a 3D vector (Å)
         *
         * Arguments:
         *      vec (list(float)):
         ******************************************************************/

        for (auto &coord : coordinates){
            coord[0] += vec[0];
            coord[1] += vec[1];
            coord[2] += vec[2];
        }
    }

    vector<unsigned long> Species::masked_atom_idxs(){
        // Atom indexes of the masked atoms in this molecule

        vector<unsigned long> indexes;

        for (unsigned long idx=0; idx<n_atoms(); idx++){
            if (atoms[idx].masked) indexes.push_back(idx);
        }

        return indexes;
    }

    array<double, 3> Species::n_vector(unsigned long i, unsigned long j){
        /*******************************************************************
         * Normalised vector between to atoms i->j
         *
         * Arguments:
         *      i (int): Atom index
         *
         *      j (int): Atom index
         *
         ******************************************************************/
        double l = distance(i, j);  // Checks for range
        array<double, 3> vec = {(coordinates[j][0] - coordinates[i][0]) / l,
                                (coordinates[j][1] - coordinates[i][1]) / l,
                                (coordinates[j][2] - coordinates[i][2]) / l};
        return vec;
    }

    void Species::print_xyz_file(const string& filename){
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

            for (int i=0; i<n_atoms(); i++){

                if (atoms[i].masked) continue;  // Skip masked atoms

                xyz_file << atoms[i].symbol   << "    "
                << coordinates[i].x() << "    "
                << coordinates[i].y() << "    "
                << coordinates[i].z() << "    " <<  '\n';
            }
            // ---------------------------------------------------

            xyz_file.close();
        }

        else throw runtime_error("Cannot open "+filename);
    }

}

