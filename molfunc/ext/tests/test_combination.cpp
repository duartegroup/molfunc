#include <iostream>
#include <fstream>
#include "molecules.h"
#include "fragments.h"
#include "combined.h"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


void print_core_xyz(){
    // xyz file appropriate for a core molecule (i.e. a dummy atom that
    // is monovalent)

    ofstream xyz_file ("core.xyz");
    if (xyz_file.is_open()){
        xyz_file << "5\n"
                    "\n"
                    "C          1.57959       -1.40470        0.00000\n"
                    "R          2.68899       -1.40471        0.00000\n"
                    "H          1.20979       -0.63118       -0.70404\n"
                    "H          1.20978       -1.18174        1.02191\n"
                    "H          1.20978       -2.40119       -0.31787\n";
        xyz_file.close();
    }
    else throw runtime_error("Unable to open core.xyz");
}


CoreMolecule core_mol(){
    print_core_xyz();
    CoreMolecule core = CoreMolecule("core.xyz");
    remove("core.xyz");

    return core;
}


CoreMolecule core_mol_two_sites(){

    ofstream xyz_file ("core.xyz");
    if (xyz_file.is_open()){
        xyz_file << "5\n"
                    "\n"
                    "C          1.57959       -1.40470        0.00000\n"
                    "R          2.68899       -1.40471        0.00000\n"
                    "R          1.20979       -0.63118       -0.70404\n"
                    "H          1.20978       -1.18174        1.02191\n"
                    "H          1.20978       -2.40119       -0.31787\n";
        xyz_file.close();
    }
    else throw runtime_error("Unable to open core.xyz");

    CoreMolecule core = CoreMolecule("core.xyz");
    remove("core.xyz");

    return core;
}


TEST_CASE("Test CombinedMolecule init from only a core"){

    auto core = core_mol();
    vector<Fragment> fragments = {};

    REQUIRE_NOTHROW(CombinedMolecule(core, fragments));
}


TEST_CASE("Test throws on unequal fragments and dummy core atoms"){

    auto core = core_mol();
    vector<Fragment> fragments = {FragmentLib::instance().fragment("Br"),
                                  FragmentLib::instance().fragment("Br")};

    // One dummy atom in the core but two fragments
    REQUIRE_THROWS(CombinedMolecule(core, fragments));
}


TEST_CASE("Test simple H3CBr combined construction") {

    auto core = core_mol();
    vector<Fragment> fragments = {FragmentLib::instance().fragment("Br")};

    auto mol = CombinedMolecule(core, fragments).to_molecule();
    int br_idx = 4;

    // Ensure the C-Br distance is reasonable
    REQUIRE((mol.distance(0, br_idx) > 1.5 && mol.distance(0, br_idx) < 2.5));

    // and that there are no short Br-H contacts
    auto h_atom_idxs = vector<int>{1, 2, 3};
    for (auto idx: h_atom_idxs){
        REQUIRE((mol.distance(idx, br_idx) > 2.0));   // r(Br-H) > 2.0 Å
    }
}


TEST_CASE("Test simple repulsive energy"){

    auto mol = CombinedMolecule(core_mol(), {});
    auto fragment = FragmentLib::instance().fragment("Br");

    // Place the fragment in a specific location
    fragment.coordinates[0] = {3.539590, -1.404700,	-0.000018};
    mol.fragments = {fragment};

    // Built molecule should have a lower repulsion than a close translation
    // of the fragment
    double rep_e = mol.repulsive_energy();

    mol.fragments[0].translate({-0.1, 0.0, 0.0});

    REQUIRE(rep_e < mol.repulsive_energy());
    rep_e = mol.repulsive_energy();

    // and if it's translated even closer
    mol.fragments[0].translate({-0.3, 0.0, 0.0});
    REQUIRE(rep_e < mol.repulsive_energy());

}


TEST_CASE("Test simple ethane combined construction") {

    auto core = core_mol();
    vector<Fragment> fragments = {FragmentLib::instance().fragment("Me")};

    auto mol = CombinedMolecule(core, fragments);
    REQUIRE(mol.repulsive_energy() < 1);
}



TEST_CASE("Test simple propane combined construction") {

    auto core = core_mol_two_sites();
    vector<Fragment> fragments = {FragmentLib::instance().fragment("Me"),
                                  FragmentLib::instance().fragment("Me")};

    auto mol = CombinedMolecule(core, fragments);

    cerr << mol.repulsive_energy() << endl;
    mol.to_molecule().print_xyz_file("tmp.xyz");

    REQUIRE(mol.repulsive_energy() < 2);
}
