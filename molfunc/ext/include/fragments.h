#ifndef MOLFUNC_EXT_FRAGMENTS_H
#define MOLFUNC_EXT_FRAGMENTS_H
#include <vector>
#include "molecules.h"


using namespace std;

namespace molfunc{

    class Fragment: public Molecule{

        public:
            vector<string> aliases;

            Fragment();
            Fragment(const string &xyz_filename);

    };


    class FragmentLib{
        // Singleton implementation from:
        // (https://stackoverflow.com/questions/86582/singleton-how-should-it-be-used)

        private:
            constexpr FragmentLib();
            // Stop the compiler generating methods of copy the object
            FragmentLib(FragmentLib const& copy);            // Not Implemented
            FragmentLib& operator=(FragmentLib const& copy); // Not Implemented

        public:
            static vector<Fragment> fragments;

            static FragmentLib& getInstance(){
                // The only instance
                // Guaranteed to be lazy initialized and destroyed correctly
                static FragmentLib instance;
                return instance;
            }
    };


}

#endif //MOLFUNC_EXT_FRAGMENTS_H
