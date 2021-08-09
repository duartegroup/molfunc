#include "atoms.h"
#include "molecules.h"
#include "stdexcept"
#include "vector"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;

bool is_close(double a, double b){
    return abs(a - b) < 1E-8;
}

bool is_close(vector<double> a, vector<double> b){

    if (a.size() != b.size()){
        throw runtime_error("Size of the two vectors not identical, "
                            "therefore not the same!");
    }

    // Ensure all elements are close
    for (int i=0; i<a.size(); i++){
        if (!is_close(a[i], b[i])){
            return false;
        }
    }

    return true;
}


TEST_CASE("Test atom must be translated by 3D vector"){

    Atom atom = Atom("H", 0.0, 0.0, 0.0);
    vector<double> xy_vec = vector<double>(2, 1.0);

    REQUIRE_THROWS(
            // Cannot translate by a 2-component vector
            atom.translate(xy_vec)
                               );

    vector<double> xyz_vec = {1.0, 0.0, 0.0};

    atom.translate(xyz_vec);
    vector<double> expected_coord = vector<double>{1.0, 0.0, 0.0};

    REQUIRE(is_close(atom.coord, expected_coord));
}


