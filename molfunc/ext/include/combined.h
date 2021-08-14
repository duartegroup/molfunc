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
    };

}



#endif //MOLFUNC_EXT_COMBINED_H
