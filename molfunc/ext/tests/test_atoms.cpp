#include "atoms.h"
#include "utils.h"
#include "stdexcept"
#include "vector"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


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

    REQUIRE(utils::is_close(atom.coord, expected_coord));
}


TEST_CASE("Test atomic numbers"){

    Atom h = Atom("H", 0.0, 0.0, 0.0);
    REQUIRE(h.atomic_number() == 1);

    Atom c = Atom("C", 0.0, 0.0, 0.0);
    REQUIRE(c.atomic_number() == 6);

    Atom og = Atom("Og", 0.0, 0.0, 0.0);
    REQUIRE(og.atomic_number() == 118);

    Atom unknown_atom = Atom("X", 0.0, 0.0, 0.0);
    REQUIRE_THROWS(unknown_atom.atomic_number());
}
