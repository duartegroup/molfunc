#include "fragments.h"
#include "utils.h"
#include "stdexcept"
#include <iostream>
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


TEST_CASE("Test fragment initialisation with aliases"){

    Fragment br = br_fragment();
    REQUIRE(br.n_atoms() == 2);

    // Should have three aliases for a bromide fragment
    REQUIRE(br.aliases.size() == 3);
}


TEST_CASE("Test fragment library construction"){

    REQUIRE(FragmentLib::instance().fragments.size() > 0);

}
