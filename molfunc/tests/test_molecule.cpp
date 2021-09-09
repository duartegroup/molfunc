#include "species/molecules.h"
#include "utils.h"
#include "vector3d.h"
#include "stdexcept"
#include "iostream"
#include <fstream>
#include "cstdio"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


void print_methane_xyz(){
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
}


Molecule methane_molecule(){
    // Generate a reasonable methane molecule

    print_methane_xyz();
    Molecule methane = Molecule("methane.xyz");
    remove("methane.xyz");

    return methane;
}


TEST_CASE("Test a molecule can be constructed from a xyz file"){

    Molecule methane = methane_molecule();

    //for (auto atom : methane.atoms){
    //    cout << atom.x() << '\t' << atom.y()  << '\t'<< atom.z() << endl;
    //}

    REQUIRE(methane.n_atoms() == 5);
    REQUIRE(utils::is_close(methane.distance(0, 1),
                            1.1,
                            0.1));

    REQUIRE(utils::is_close(methane.coordinates[0][0], 1.57959));

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

    vector<Atom3D> atoms = {Atom3D("Os", 0.0, 0.0, 0.0),
                            Atom3D("Sn", 2.5, 0.0, 0.0)};
    Molecule mol = Molecule(atoms);

    // The twp atoms are bonded
    REQUIRE(mol.graph.n_neighbours(0) == 1);

}


TEST_CASE("Test throws on printing an xyz file with no atoms"){
    Molecule mol = Molecule();

    REQUIRE_THROWS(mol.print_xyz_file("tmp.xyz"));
}


TEST_CASE("Test core molecule construction"){

    print_methane_xyz();
    CoreMolecule mol = CoreMolecule("methane.xyz",
                                    {1});
    REQUIRE(mol.n_atoms() == 5);

    // Deleting atom 1 (a hydrogen) will mask it
    REQUIRE(mol.n_unmasked_atoms() == 4);
    REQUIRE(mol.masked_atom_idxs().size() == 1);
    REQUIRE(mol.masked_atom_idxs()[0] == 1);

    // Cannot construct a core molecule where the deleted atom(s) is(are) not
    // monovalent
    REQUIRE_THROWS(CoreMolecule("methane.xyz", {0}));

    // nor if the deleted atom is out of range
    REQUIRE_THROWS(CoreMolecule("methane.xyz", {5}));

    // nor if the .xyz file doesn't contain any dummy atoms and no atoms_to_del
    // are specified
    REQUIRE_THROWS(CoreMolecule("methane.xyz"));

    remove("methane.xyz");
}


TEST_CASE("Test molecule translation"){
    Molecule methane = methane_molecule();

    Vector3D vec = {-methane.coordinates[0][0],
                    -methane.coordinates[0][1],
                    -methane.coordinates[0][2]};
    methane.translate(vec);
    methane.print_xyz_file("tmp.xyz");

    Molecule regen_methane = Molecule("tmp.xyz");

    // Translating by the negative of the Origin-C vecotr should
    // leave the carbon at (0, 0, 0) i.e. the origin
    REQUIRE(utils::is_close(regen_methane.coordinates[0], { 0, 0, 0 }));

    remove("tmp.xyz");
}


TEST_CASE("Test coordinate vectors"){

    auto mol = methane_molecule();

    auto vec = mol.coordinates[0] - mol.coordinates[1];

    REQUIRE(utils::is_close(vec.length(),
                            1.1,
                            0.1));

    // Normalising should make it very close to a unit vector
    vec.normalise();
    REQUIRE(utils::is_close(vec.length(), 1.0));
}
