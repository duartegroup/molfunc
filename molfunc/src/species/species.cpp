#include "species/species.h"
#include "fstream"
#include "string"
#include "iomanip"
#include "cmath"
#include "memory"
#include "iostream"


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

        if (i >= n_atoms() || j >= n_atoms() || coordinates.empty()){
            throw runtime_error("Invalid index: "+to_string(i)+" or "+ to_string(j)+
                                " must be present in the molecule");
        }

        double sq_dist = 0.0;

        for (auto k=0; k<3; k++){
            double tmp = coordinates[i][k] - coordinates[j][k];
            sq_dist += tmp * tmp;
        }

        return sqrt(sq_dist);
    }

    void Species::translate(const Vector3D &vec){
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

    Vector3D Species::n_vector(unsigned long i, unsigned long j){
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

        return (coordinates[j] - coordinates[i]) / l;
    }

    void Species::print_xyz_file(const string& filename,
                                 bool append){
        /*********************************************************
         * Generate a standard .xyz file for a molecule
         *
         * Arguments:
         *      filename (str):
         ********************************************************/

        if (atoms.empty()){
            throw runtime_error("Could not print a .xyz file- had no atoms");
        }

        ofstream xyz_file;
        if (append) xyz_file.open(filename, ios_base::app);
        else xyz_file.open(filename);

        if (xyz_file.is_open()){

            xyz_file << fixed;
            xyz_file << setprecision(6);

            // ---------------------------------------------------
            xyz_file << to_string(n_unmasked_atoms()) << '\n'
            << "molfunc generated" << '\n';

            for (unsigned long i=0; i<n_atoms(); i++){

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

    void Species::print_xyz_file(const string& filename){
        // Print the .xyz file of this molecule, will override
        print_xyz_file(filename, false);
    }

    void Species::append_xyz_file(const string &filename) {
        // Append to a (possibly) existing .xyz file
        print_xyz_file(filename, true);
    }

    void Species::rotate(const RotationMatrix &R) {
        /*********************************************************
         * Rotate this species using a rotation matrix
         *
         * Arguments:
         *      R (RotationMatrix): (3x3) Matrix that applies the
         *                          rotation
         ********************************************************/

        for (auto &coord : coordinates){
            double x = coord[0];
            double y = coord[1];
            double z = coord[2];

            coord[0] = R[0][0] * x + R[0][1] * y + R[0][2] * z;
            coord[1] = R[1][0] * x + R[1][1] * y + R[1][2] * z;
            coord[2] = R[2][0] * x + R[2][1] * y + R[2][2] * z;
        }
    }

    void Species::rotate(const RotationMatrix &R, const Coordinate &origin){
        /*********************************************************
         * Rotate this species using a rotation matrix about a
         * defined origin in space
         *
         * Arguments:
         *      R (RotationMatrix): (3x3)
         *
         *      origin (Coordinate):
         ********************************************************/

        for (auto &coord : coordinates) coord -= origin;
        rotate(R);
        for (auto &coord : coordinates) coord += origin;
    }

    void Species::rotate(const RotationMatrix &R, unsigned long atom_idx){
        rotate(R, Coordinate(coordinates[atom_idx]));
    }

    unsigned long Species::no_masked_idx(unsigned long idx) {
        /*************************************************
         * Convert an atom index where the all atoms
         * are present to one where the masked/dummy
         * atoms have been deleted e.g. for atoms
         *
         *      [[C, 0.0, 0.0, 0.0],
         *      [R, 0.0, 0.0, 0.0],
         *      [C, 0.0, 0.0, 0.0]]
         *
         *  then: no_masked_idx(0) -> 0
         *        no_masked_idx(2) -> 1
         *
         ************************************************/
        unsigned long n_masked_atoms = 0;

        for (unsigned long i=0; i<n_atoms(); i++){
            if (atoms[i].masked){
                if (i == idx){
                    throw runtime_error("Cannot index a "
                                        "masked atom");
                }
                n_masked_atoms += 1;
                continue;
            }

            if (i == idx) return idx - n_masked_atoms;
        }

        throw out_of_range("No valid index present");
    }

}

