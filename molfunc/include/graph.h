#ifndef MOLFUNC_EXT_GRAPH_H
#define MOLFUNC_EXT_GRAPH_H
#include "iostream"
#include "string"
#include "vector"
#include <unordered_map>


using namespace std;

namespace molfunc{

    class Node{

        public:
            unsigned long idx = 0;
            string symbol = "X";
            vector<unsigned long> neighbours;

            Node();
            Node(unsigned long u,
                 string symbol);

            bool operator < (const Node& node) const{
                // Less than operator for comparison of two nodes
                return (idx < node.idx);
            }
    };

    class Edge{

        public:
            unsigned long u = 0;
            unsigned long v = 1;

            Edge();
            Edge(unsigned long u,
                 unsigned long v);
    };

    class Graph {

        protected:
            unordered_map<unsigned long, Node> nodes;
            vector<Edge> edges;

        public:

            explicit Graph();

            void add_node(const Node &node);
            unsigned long n_nodes();

            void add_edge(unsigned long u,
                          unsigned long v);
            void add_edge(const Edge &edge);
            unsigned long n_edges();

            unsigned long first_neighbour(unsigned long u);
            unsigned long first_non_dummy_neighbour(unsigned long u);
            unsigned long n_neighbours(unsigned long u);

            string atomic_symbol(unsigned long u);
    };

}


#endif //MOLFUNC_EXT_GRAPH_H
