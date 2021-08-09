#include <fstream>
#include "string"
#include "molecules.h"
#include "utils.h"


using namespace std;

namespace molfunc{

    Molecule::Molecule() = default;

    Molecule::Molecule(string xyz_filename) {
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
         *      // mol.n_atoms() = 5;
         ********************************************************/

        string line, item;
        ifstream xyz_file (xyz_filename);

        int decl_n_atoms = 0;

        // Iterate through the xyz file
        while (getline(xyz_file, line, '\n')) {

            if (line.empty()) {   // Ignore any blank lines etc.
                continue;
            }

            // Assign the number of declared atoms
            if (decl_n_atoms == 0) {
                decl_n_atoms = stoi(line);
                continue;
            }

            vector<string> xyz_items = utils::split(line, ' ');
            if (xyz_items.size() != 4){
                throw runtime_error("Malformatted xyz file, expecting a *A  0.0 0.0 0.0* structure");
            }

            atoms.emplace_back(xyz_items[0],         // Atomic symbol
                               stod(xyz_items[1]),   // x
                               stod(xyz_items[2]),   // y
                               stod(xyz_items[3]));  // z

        }
        xyz_file.close();

        if (decl_n_atoms != n_atoms()) {
            throw runtime_error("Number of atoms "+to_string(n_atoms())+" not equal to the number declared "
                                + to_string(decl_n_atoms));
        }
    }


    unsigned long Molecule::n_atoms() {
        return atoms.size();
    }

}

