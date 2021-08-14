#ifndef MOLFUNC_EXT_FRAGMENTS_H
#define MOLFUNC_EXT_FRAGMENTS_H
#include <vector>
#include <string>
#include "atoms.h"
#include "molecules.h"


using namespace std;

namespace molfunc{

    class Fragment: public Molecule{

        public:
            vector<string> aliases;

            // Constructors
            Fragment();
            explicit Fragment(const string& xyz_filename);
            Fragment(const vector<Atom>& atoms,
                     const string& title);
    };


    class FragmentLib{
        // Singleton implementation from:
        // (https://stackoverflow.com/questions/86582/singleton-how-should-it-be-used)

        private:
            FragmentLib();

        public:
            FragmentLib(FragmentLib const& copy) = delete;            // Not Implemented
            FragmentLib& operator=(FragmentLib const& copy) = delete; // Not Implemented

            vector<Fragment> fragments;

            static FragmentLib& instance(){
                // The only instance
                // Guaranteed to be lazy initialized and destroyed correctly
                static FragmentLib instance;
                return instance;
            }
    };


}

#endif //MOLFUNC_EXT_FRAGMENTS_H