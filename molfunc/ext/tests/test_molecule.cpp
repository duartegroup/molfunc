#include "atoms.h"
#include "molecules.h"
#include "stdexcept"
#include "vector"
#include <iostream>
#include <fstream>
#include "cstdio"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


TEST_CASE("Test a molecule can be constructed from a xyz file"){

    ofstream xyz_file ("methane.xyz");
    if (xyz_file.is_open()){
        xyz_file << "5\n"
                    "\n"
                    "C          1.57959       -1.40470        0.00000\n"
                    "H          2.68899       -1.40471        0.00000\n"
                    "H          1.20979       -0.63118       -0.70404\n"
                    "H          1.20978       -1.18174        1.02191\n"
                    "H          1.20978       -2.40119       -0.31787\n";
        xyz_file.close();
    }
    else {
        cout << "Unable to open methane.xyz";
        return;
    }

    Molecule methane = Molecule("methane.xyz");

    REQUIRE(methane.n_atoms() == 5);
    remove("methane.xyz");
}




