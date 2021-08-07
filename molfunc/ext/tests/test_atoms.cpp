#include "atoms.h"
#include "vector"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;

TEST_CASE("Test atom must be translated by 3D vector"){

    Atom atom = Atom("H", 0.0, 0.0, 0.0);
    REQUIRE_THROWS(
            atom.translate(vector<double>(2, 1.0))
                               );
}


