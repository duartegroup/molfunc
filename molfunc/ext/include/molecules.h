#ifndef MOLFUNC_MOLECULES_H
#define MOLFUNC_MOLECULES_H
#include "atoms.h"
#include "species.h"
#include "graph.h"
#include "array"
#include "vector"
#include "string"


using namespace std;


namespace molfunc{

    class Molecule: public Species{

        public:

            molfunc::Graph graph;

            Molecule();
            explicit Molecule(const string& xyz_filename);
            explicit Molecule(const vector<Atom>& atoms);

            void construct_graph();

            void assign_coordinates();

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
