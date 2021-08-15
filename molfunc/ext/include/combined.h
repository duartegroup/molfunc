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
            CombinedMolecule(CoreMolecule &core,
                             vector<Fragment> &fragments);

            void build();

            Molecule to_molecule();


        protected:
            void translate_fragment(Fragment &fragment,
                                    unsigned long dummy_atom_idx);
    };

}



#endif //MOLFUNC_EXT_COMBINED_H
