#include "atoms.h"
#include <utility>
#include <algorithm>
#include <vector>
#include <stdexcept>

using namespace std;


namespace molfunc{

    // Default constructors are required for Cython wrapping
    Atom::Atom() = default;

    Atom::Atom(string symbol, double x, double y, double z) {
        /*********************************************************
         * Construct an atom from an atomic symbol and its position
         * in 3D space
         *
         * Arguments:
         *      symbol (string): Atomic symbol
         *      x (float):       Cartesian x coordinate (Å)
         *      y (float):                 y coordinate (Å)
         *      z (float):                 z coordinate (Å)
         ********************************************************/

        this->symbol = std::move(symbol);
        this->coord = vector<double>{x, y, z};

    }

    void Atom::translate(vector<double> &vec){
        /*********************************************************
         * Translate an atom by a 3-component vector (Å)
         *
         * Arguments:
         *      vec (vector): Shift vector in 3D space
         *
         *  Raises:
         *      (runtime_error):
         ********************************************************/
        if (vec.size() != 3){
            throw runtime_error("Must have vector of length 3 to translate");
        }

        for (int i=0; i<3; i++){
            coord[i] += vec[i];
        }
    }

    unsigned int Atom::atomic_number(){
        /*********************************************************
         * Determine the atomic number of this element by its
         * symbol.
         *
         * Returns:
         *      (int): Atomic number [1, 118]
         *
         * Raises:
         *      (runtime_error): If the symbol cannot be found in
         *                    the periodic table of known elements
         ********************************************************/

        // Defines a special dummy atom ("R") with atomic number 0
        vector<string> elements = {"R", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
                                   "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
                                   "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
                                   "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                                   "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
                                   "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
                                   "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
                                   "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                                   "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
                                   "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
                                   "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

        for (unsigned int i=0; i<elements.size(); i++){
            if (elements[i] == symbol){
                return i;
            }
        }

        throw runtime_error("Failed to find atomic number for "+symbol);
    }

    bool Atom::is_dummy(){return (atomic_number() == 0);}

}
