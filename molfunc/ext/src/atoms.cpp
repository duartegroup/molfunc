#include "atoms.h"
#include <utility>
#include <algorithm>
#include <array>
#include <vector>
#include <stdexcept>
#include "iostream"

using namespace std;


namespace molfunc{

    // Default constructors are required for Cython wrapping
    Atom::Atom() = default;

    Atom::Atom(const string& symbol, double x, double y, double z) {
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

        this->symbol = symbol;
        this->coord = {x, y, z};

        if (symbol == "R") masked = true;

    }

    void Atom::translate(array<double, 3> &vec){
        /*********************************************************
         * Translate an atom by a 3-component vector (Å)
         *
         * Arguments:
         *      vec (vector): Shift vector in 3D space
         *
         *  Raises:
         *      (runtime_error):
         ********************************************************/

        for (int i=0; i<3; i++){
            coord[i] += vec[i];
        }
    }

    unsigned int Atom::atomic_number() const{
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

    bool Atom::is_dummy() const{return (symbol == "R");}

    double Atom::x(){
        if (ptr_x != nullptr) coord[0] = *ptr_x;
        return coord[0];
    }

    double Atom::y(){
        if (ptr_y != nullptr) coord[1] = *ptr_y;
        return coord[1];
    }

    double Atom::z(){
        if (ptr_z != nullptr) coord[2] = *ptr_z;
        return coord[2];
    }
}
