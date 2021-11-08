#ifndef MOLFUNC_EXT_FRAGMENTS_H
#define MOLFUNC_EXT_FRAGMENTS_H
#include "vector"
#include "string"
#include "iostream"
#include "atoms.h"
#include "grid.h"
#include "molecules.h"


using namespace std;

namespace molfunc{

    class Fragment: public Molecule{

        protected:
            vector<Coordinate> cached_coordinates;
            RotationMatrix rotation_matrix = RotationMatrix();

        public:
            vector<string> aliases;
            unsigned long dummy_idx = 0;
            unsigned long dummy_nn_idx = 0;
            Grid3D rot_grid_w = Grid3D(0.001, 3.14, 10);

            // Constructors
            Fragment();
            explicit Fragment(const string& xyz_filename);
            Fragment(const Fragment &fragment);
            Fragment(vector<Atom3D> atoms, vector<string> aliases);


            void reset_coordinates();
            void cache_coordinates();

            void rotate(GridPoint &grid_point);
            void rotate_about_dummy_nn(GridPoint &grid_point);
    };


    class FragmentLib{   // Singleton

        private:
            FragmentLib();

        public:
            FragmentLib(FragmentLib const& copy) = delete;            // Not Implemented
            FragmentLib& operator=(FragmentLib const& copy) = delete; // Not Implemented

            static FragmentLib& instance(){
                // The only instance
                // Guaranteed to be lazy initialized and destroyed correctly
                static FragmentLib instance;
                return instance;
            }

            vector<Fragment> fragments;

            Fragment fragment(const string& name);

            vector<vector<Fragment>> fragments_n_repeats(unsigned long n);
    };


}

#endif //MOLFUNC_EXT_FRAGMENTS_H
