#ifndef MOLFUNC_MOLECULES_H
#define MOLFUNC_MOLECULES_H
#include "atoms.h"
#include "graph.h"
#include "array"
#include "vector"
#include "string"

using namespace std;

namespace molfunc{

    class Molecule {

        public:

            molfunc::Graph graph;

            vector<array<double, 3>> coordinates;

            Molecule();
            explicit Molecule(const string& xyz_filename);
            explicit Molecule(const vector<Atom>& atoms);

            unsigned long n_atoms() const;
            unsigned long n_masked_atoms();
            unsigned long n_unmasked_atoms();

            double distance(unsigned long i, unsigned long j);

            void construct_graph();

            void assign_coordinates();

            void print_xyz_file(const string& filename);

        protected:
            string xyz_title_line;
            vector<Atom> atoms;

        void set_atoms(const string& xyz_filename);
            bool is_bonded_on_distance(unsigned long i, unsigned long j);

    };


    class CoreMolecule: public Molecule{

        public:
            CoreMolecule(const string& xyz_filename,
                         const vector<unsigned int>& atoms_to_del);

            CoreMolecule(const string& xyz_filename);

    };
}

#endif //MOLFUNC_MOLECULES_H
