#ifndef MOLFUNC_EXT_SPECIES_H
#define MOLFUNC_EXT_SPECIES_H
#include "vector"
#include "atoms.h"


namespace molfunc{

    class Species {

        public:
            Species();

            vector<Coordinate> coordinates;
            vector<Atom> atoms;

            unsigned long n_atoms() const;
            unsigned long n_masked_atoms();
            unsigned long n_unmasked_atoms();

            double distance(unsigned long i, unsigned long j);

            void translate(array<double, 3> vec);

            array<double, 3> n_vector(unsigned long i, unsigned long j);

            vector<unsigned long> masked_atom_idxs();

            void print_xyz_file(const string& filename);

    };
}


#endif //MOLFUNC_EXT_SPECIES_H
