#include "utils.h"
#include "grid.h"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


TEST_CASE("Test simple -1 -> 1 grid"){
    auto grid = Grid3D(-1, 1, 3);

    // First element should be is (-1, -1, -1)
    for (const auto &p : grid[0]){
        REQUIRE(utils::is_close(p, -1));
    }

    // and the final element is (1, 1, 1)
    for (const auto &p : grid.back()){
        REQUIRE(utils::is_close(p, 1));
    }
}

