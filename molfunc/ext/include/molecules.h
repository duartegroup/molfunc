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
            Molecule(string xyz_filename);

            unsigned long n_atoms();

            void print_xyz_file();

    };
}

#endif //MOLFUNC_MOLECULES_H
