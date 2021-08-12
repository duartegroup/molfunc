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
            explicit Molecule(const string& xyz_filename);

            unsigned long n_atoms();
            double distance(unsigned long i, unsigned long j);

            void construct_graph();

            void print_xyz_file(const string& filename);

        protected:
            string xyz_title_line;

            void set_atoms(const string& xyz_filename);
            bool is_bonded_on_distance(unsigned long i, unsigned long j);

    };
}

#endif //MOLFUNC_MOLECULES_H
