#include "atoms.h"
#include <utility>
#include "vector"
#include "stdexcept"

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


}
