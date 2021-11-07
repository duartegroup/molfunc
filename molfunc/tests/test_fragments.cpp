#include "species/fragments.h"
#include "rotation.h"
#include "utils.h"
#include "stdexcept"
#include "iostream"
#include <fstream>
#include "cstdio"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


Fragment br_fragment(){

    ofstream xyz_file ("br.xyz");
    if (xyz_file.is_open()){
        xyz_file << "2\n"
                    "Br[*] bromide,bromo,br\n"
                    "Br  0.92450   0.00000   0.00000\n"
                    "R   -0.92450  0.00000   0.00000\n";
        xyz_file.close();
    }
    else {
        throw runtime_error("Unable to open br.xyz");
    }
    Fragment br = Fragment("br.xyz");
    remove("br.xyz");

    return br;
}


void print_br_xyz(){
    ofstream xyz_file ("br.xyz");
    if (xyz_file.is_open()){
        xyz_file << "1\n"
                    "\n"
                    "Br          0.0       0.0        0.0\n";
        xyz_file.close();
    }
    else throw runtime_error("Unable to open br.xyz");
}


TEST_CASE("Test fragment initialisation with aliases"){

    Fragment br = br_fragment();
    REQUIRE(br.n_atoms() == 2);

    // Should have three aliases for a bromide fragment
    REQUIRE(br.aliases.size() == 3);
}


TEST_CASE("Test fragment library construction"){

    REQUIRE_FALSE(FragmentLib::instance().fragments.empty());

}


TEST_CASE("Test retrieval from fragment library"){

    auto frag = FragmentLib::instance().fragment("Br");

    REQUIRE(frag.n_atoms() == 2);  // Br and R atom
    REQUIRE(frag.n_unmasked_atoms() == 1); // Br

    REQUIRE_THROWS(FragmentLib::instance().fragment("Non-existing-frag"));
}


TEST_CASE("Test fragment constructor throws with no dummy atoms"){

    print_br_xyz();
    REQUIRE_THROWS(Fragment("br.xyz"));
    remove("br.xyz");
}


TEST_CASE("Test multiple fragments from lib"){

    Fragment br_1 = FragmentLib::instance().fragment("Br");
    auto coord = br_1.coordinates[0];

    br_1.translate({1.0, 0.0, 0.0});

    // Should hava translated the atom along the x axis
    REQUIRE(!utils::is_close(coord[0], br_1.coordinates[0][0]));

    // but a new Br fragment should not have also been translated
    Fragment br_2 = FragmentLib::instance().fragment("Br");

    REQUIRE(utils::is_close(coord[0], br_2.coordinates[0][0]));
}


TEST_CASE("Test coordinate reset"){
    auto fragment = br_fragment();
    fragment.cache_coordinates();

    // Apply a random ish rotation to the fragment
    GridPoint point = {1.0, 0.1, 0.2};
    fragment.rotate(point);

    REQUIRE_FALSE(utils::is_close(fragment.coordinates[0].x(), 0.92450));

    // Resetting the coordinates should return them to the original values
    fragment.reset_coordinates();
    REQUIRE(utils::is_close(fragment.coordinates[0].x(), 0.92450));
}


TEST_CASE("Test repeated fragment combinations"){
    /* Should be able to generate all possible
     of fragments with a defined number of repeats
     e.g a fragment library of size 2 should generate:
     n = 1   --> [[frag1, frag2]]
     n = 2   --> [[frag1, frag1],
                  [frag1, frag2],
                  [frag2, frag1],
                  [frag2, frag2]]
    */

    auto n_fragments = FragmentLib::instance().fragments.size();

    auto vec_fragments = FragmentLib::instance().fragments_n_repeats(1);
    REQUIRE(vec_fragments.size() == n_fragments);
    REQUIRE(vec_fragments[0].size() == 1);

    vec_fragments = FragmentLib::instance().fragments_n_repeats(2);

    REQUIRE(vec_fragments.size() == n_fragments*n_fragments);

    // Don't generate > 24 million combinations..
    REQUIRE_THROWS(FragmentLib::instance().fragments_n_repeats(5));
}

