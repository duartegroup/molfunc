#include "graph.h"
#include "stdexcept"
#include "iostream"
#include <utility>

using namespace std;

namespace molfunc{

    Node::Node() = default;

    Node::Node(unsigned long u, string symbol) {
        /*********************************************************
         * Construct a node of a molecular graph using
         * an index, and the atomic symbol
         *
         * Arguments:
         *      idx (int): Atom index e.g. 0
         *
         *      symbol (str): Atomic symbol e.g. "H"
         *
         * Example:
         *      Node node = Node(0, "H")
         ********************************************************/

        this->idx = u;
        this->symbol = std::move(symbol);
    }

    Edge::Edge() = default;

    Edge::Edge(unsigned long u, unsigned long v) {
        /*********************************************************
         * Edge between two nodes in a graph, idx and v
         *
         * Arguments:
         *      idx (int): Atom index
         *
         *      v (int): Atom index
         *
         * Example:
         *      Edge edge = Edge(0, 1)
         ********************************************************/
        this->u = u;
        this->v = v;
    }

    Graph::Graph() = default;

    void Graph::add_node(const Node &node) {
        /*********************************************************
         * Add a node to this molecular graph if it is already
         * present then throw a warning but do not append it again
         *
         * Arguments:
         *      idx (int): Index of the node
         ********************************************************/
        nodes[node.idx] = node;
    }

    unsigned long Graph::n_nodes() {
        // Number of nodes in the graph
        return nodes.size();
    }

    void Graph::add_edge(unsigned long u, unsigned long v) {
        /*********************************************************
         * Add an edge between two atoms in a molecule graph,
         * they must both exist
         *
         * Arguments:
         *      idx (int): Atom index
         *
         *      v (int): Atom index
         *
         *  Raises:
         *      (runtime_error): If the two atom indexes are
         *                       not present in the graph, or
         ********************************************************/

        if (n_nodes() < 2){
            throw runtime_error("Cannot add an edge for a graph with "
                                "fewer than two nodes");
        }

        if (u == v){
            throw runtime_error("Self referring are not supported");
        }

        if (nodes.find(u) != nodes.end()){
            throw runtime_error("Node "+to_string(u)+" not present, cannot "
                                "add edge");
        }

        if (nodes.find(v) != nodes.end()){
            throw runtime_error("Node "+to_string(v)+" not present, cannot "
                                                     "add edge");
        }

        edges.emplace_back(u, v);

        nodes[u].neighbours.push_back(v);
        nodes[v].neighbours.push_back(u);
    }

    void Graph::add_edge(const Edge &edge) {
        // Add an edge
        return add_edge(edge.u, edge.v);
    }

    unsigned long Graph::n_edges() {
        // Number of edges in the graph
        return edges.size();
    }

}
