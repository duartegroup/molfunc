#include "molecules.h"
#include "utils.h"
#include "stdexcept"
#include <iostream>
#include <fstream>
#include "cstdio"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


Molecule methane_molecule(){
    // Generate a reasonable methane molecule

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
    else throw runtime_error("Unable to open methane.xyz");

    Molecule methane = Molecule("methane.xyz");
    remove("methane.xyz");

    return methane;
}


TEST_CASE("Test a molecule can be constructed from a xyz file"){

    Molecule methane = methane_molecule();

    REQUIRE(methane.n_atoms() == 5);
    REQUIRE(utils::is_close(methane.distance(0, 1),
                            1.1,
                            0.1));

    // Printing the molecule should be able to be read (i.e. be valid)
    methane.print_xyz_file("tmp.xyz");

    Molecule regen_methane = Molecule("tmp.xyz");
    REQUIRE(regen_methane.n_atoms() == 5);
    remove("tmp.xyz");
}


TEST_CASE("Test throws if xyz file does not exist"){
    REQUIRE_THROWS(Molecule("x.xyz"));
}

TEST_CASE("Test throws if not an .xyz file"){

    ofstream xyz_file ("tmp.txt");
    if (xyz_file.is_open()) {
        xyz_file << "something" << '\n';
        xyz_file.close();

        REQUIRE_THROWS(Molecule("tmp.txt"));
        remove("tmp.txt");
    }
    else throw runtime_error("Unable to open tmp.txt");

}


TEST_CASE("Test a graph is constructed for a molecule"){

    Molecule methane = methane_molecule();

    REQUIRE(methane.graph.n_nodes() == 5);
    REQUIRE(methane.graph.n_edges() == 4);

    REQUIRE(methane.graph.n_neighbours(0) == 4);
    REQUIRE(methane.graph.n_neighbours(1) == 1);

}


TEST_CASE("Test bond definitions"){

    Molecule mol = Molecule();
    mol.atoms.push_back(Atom("Os", 0.0, 0.0, 0.0));
    mol.atoms.push_back(Atom("Sn", 2.5, 0.0, 0.0));
    mol.construct_graph();

    // The twp atoms are bonded
    REQUIRE(mol.graph.n_neighbours(0) == 1);

}