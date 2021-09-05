#ifndef MOLFUNC_EXT_COMBINED_H
#define MOLFUNC_EXT_COMBINED_H
#include "molecules.h"
#include "fragments.h"


namespace molfunc{

    class CombinedMolecule {

        public:

            CoreMolecule core;
            vector<Fragment> fragments;

            CombinedMolecule();
            CombinedMolecule(CoreMolecule core,
                             vector<Fragment> fragments);

            double repulsive_energy();
            double repulsive_energy(const Fragment& fragment);
            void build();

            Molecule to_molecule();


        protected:
            void rotate_fragments_global();
            void exclude_rotational_space(Fragment &fragment,
                                          double threshold);

            void translate_fragment(Fragment &fragment,
                                    unsigned long Ra_idx);

            void centre_core_and_fragment_to(Fragment &fragment);
    };

}


#endif //MOLFUNC_EXT_COMBINED_H