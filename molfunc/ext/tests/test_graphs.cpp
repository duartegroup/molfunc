#include "graph.h"
#include "catch2/catch.hpp"

using namespace std;
using namespace molfunc;


TEST_CASE("Test a graph can be constructed from nodes and edges"){

    Graph graph = Graph();

    REQUIRE(graph.n_nodes() == 0);
    REQUIRE(graph.n_edges() == 0);

    graph.add_node(Node(0, "X"));
    REQUIRE(graph.atomic_symbol(0) == "X");

    graph.add_node(Node(1, "X"));
    graph.add_edge(0, 1);

    REQUIRE(graph.n_edges() == 1);
    REQUIRE(graph.n_neighbours(0) == 1);
    REQUIRE(graph.n_neighbours(1) == 1);

}


TEST_CASE("Test a graph addition exceptions"){

    Graph graph = Graph();
    graph.add_node(Node(0, "X"));

    // Cannot have an atom symbol for a non-existing node
    REQUIRE_THROWS(graph.atomic_symbol(1));

    // Cannot add an edge where atom 1 does not exist
    REQUIRE_THROWS(graph.add_edge(0, 1));

    // or between the same node
    REQUIRE_THROWS(graph.add_edge(0, 0));

    graph.add_node(Node(2, "X"));

    // or where there are enough edges but not the correct ones
    REQUIRE_THROWS(graph.add_edge(Edge(0, 1)));

}

