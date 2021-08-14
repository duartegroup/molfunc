#include "atoms.h"
#include "utils.h"
#include "array"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


TEST_CASE("Test atom must be translated by 3D vector"){

    Atom atom = Atom("H", 0.0, 0.0, 0.0);

    array<double, 3> xyz_vec = {1.0, 0.0, 0.0};

    atom.translate(xyz_vec);
    vector<double> expected_coord = vector<double>{1.0, 0.0, 0.0};

    REQUIRE(utils::is_close(atom.x(), expected_coord[0]));
    REQUIRE(utils::is_close(atom.y(), expected_coord[1]));
    REQUIRE(utils::is_close(atom.z(), expected_coord[2]));

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
