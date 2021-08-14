#include "species.h"
#include <fstream>
#include <string>
#include <iomanip>
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

}

