#ifndef MOLFUNC_EXT_SPECIES_H
#define MOLFUNC_EXT_SPECIES_H
#include "vector"
#include "atoms.h"


namespace molfunc{

    class Species {

        public:
            Species();

            vector<array<double, 3>> coordinates;
            vector<Atom> atoms;

            unsigned long n_atoms() const;
            unsigned long n_masked_atoms();
            unsigned long n_unmasked_atoms();

            void print_xyz_file(const string& filename);

    };
}


#endif //MOLFUNC_EXT_SPECIES_H
