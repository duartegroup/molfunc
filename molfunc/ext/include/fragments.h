#ifndef MOLFUNC_EXT_FRAGMENTS_H
#define MOLFUNC_EXT_FRAGMENTS_H
#include <vector>
#include <string>
#include "atoms.h"
#include "grid.h"
#include "molecules.h"


using namespace std;

namespace molfunc{

    class Fragment: public Molecule{

        protected:
            vector<Coordinate> cached_coordinates;

        public:
            vector<string> aliases;
            unsigned long dummy_idx = 0;
            Grid3D rot_grid_w = Grid3D(0.001, 3.14, 10);

            // Constructors
            Fragment();
            explicit Fragment(const string& xyz_filename);
            Fragment(const Fragment &fragment);

            void reset_coordinates();

    };


    class FragmentLib{   // Singleton

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

            Fragment fragment(const string& name);
    };


}

#endif //MOLFUNC_EXT_FRAGMENTS_H
