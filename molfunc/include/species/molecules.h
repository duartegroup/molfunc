#ifndef MOLFUNC_MOLECULES_H
#define MOLFUNC_MOLECULES_H
#include "array"
#include "vector"
#include "string"
#include "iostream"
#include "atoms.h"
#include "species.h"
#include "graph.h"


using namespace std;


namespace molfunc{

    class Molecule: public Species{

        public:

            molfunc::Graph graph;

            Molecule();
            explicit Molecule(const string& xyz_filename);
            explicit Molecule(vector<Atom3D>& atom);

            void construct_graph();

        protected:
            string xyz_title_line;

            void set_atoms(const string& xyz_filename);
            bool is_bonded_on_distance(unsigned long i, unsigned long j);

    };


    class CoreMolecule: public Molecule{

        public:
            CoreMolecule();
            CoreMolecule(const string& xyz_filename,
                         const vector<unsigned int>& atoms_to_del);
            explicit CoreMolecule(const string& xyz_filename);

    };
}

#endif //MOLFUNC_MOLECULES_H
