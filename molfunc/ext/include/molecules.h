#ifndef MOLFUNC_MOLECULES_H
#define MOLFUNC_MOLECULES_H
#include "atoms.h"
#include "graph.h"
#include "string"

using namespace std;

namespace molfunc{

    class Molecule {

        public:

            molfunc::Graph graph;
            vector<Atom> atoms;

            Molecule();
            Molecule(const string& xyz_filename);

            unsigned long n_atoms() const;
            double distance(unsigned long i, unsigned long j);

            void print_xyz_file();

        protected:
            void set_atoms(const string& xyz_filename);
            void construct_graph();
            bool is_bonded_on_distance(unsigned long i, unsigned long j);

    };
}

#endif //MOLFUNC_MOLECULES_H
