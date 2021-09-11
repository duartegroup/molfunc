#include "atoms.h"
#include "utils.h"
#include "array"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


TEST_CASE("Test atomic numbers"){

    auto h = Atom3D("H", 0.0, 0.0, 0.0);
    REQUIRE(h.atomic_number() == 1);

    auto c = Atom3D("C", 0.0, 0.0, 0.0);
    REQUIRE(c.atomic_number() == 6);

    auto og = Atom3D("Og", 0.0, 0.0, 0.0);
    REQUIRE(og.atomic_number() == 118);

    auto unknown_atom = Atom3D("X", 0.0, 0.0, 0.0);
    REQUIRE_THROWS(unknown_atom.atomic_number());
}
