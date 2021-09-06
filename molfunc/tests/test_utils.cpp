#include "string"
#include "utils.h"
#include "catch2/catch.hpp"


using namespace std;
using namespace molfunc;

TEST_CASE("Test can lower case a string"){

    string test = "TeSt";
    REQUIRE(utils::to_lower(test)  == "test");

}
