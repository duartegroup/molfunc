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


TEST_CASE("Test minimum energy point no grid points"){

    auto grid = Grid3D(-1.0, 1.0, 0);

    // No grid points to evaluate on
    REQUIRE_THROWS(grid.minimum_energy_point());

}


TEST_CASE("Test minimum energy point"){

    auto grid = Grid3D(-1.0, 1.0, 3);
    for (auto &point : grid){
        point.energy = 1.0;
    }
    grid[1].energy = 0.1;

    auto point = grid.minimum_energy_point();

    // Should pick the one with the smallest energy...
    REQUIRE(utils::is_close(point.energy, 0.1));
    // and be the second point in the list
    REQUIRE(utils::is_close(point[0], -1.0));
    REQUIRE(utils::is_close(point[1], -1.0));
    REQUIRE(utils::is_close(point[2], 0.0));

}

