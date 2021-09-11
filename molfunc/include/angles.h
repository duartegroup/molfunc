#ifndef MOLFUNC_EXT_ANGLES_H
#define MOLFUNC_EXT_ANGLES_H
#include "vector"
#include "iostream"
#include "stdexcept"
#include "atoms.h"


namespace molfunc{

    class AnglePotential {

        public:
            double cos_phi0 = 0.0;
            double half_k = 1.0;

            array<unsigned long, 3> atom_idxs = {0, 1, 2};

            AnglePotential();
            AnglePotential(double phi0, double k);
            AnglePotential(unsigned long idx0,
                           unsigned long idx1,
                           unsigned long idx2,
                           double phi0, double k);

            double value(vector<Coordinate> &coordinates);

    };


    class AnglePotentials: public vector<AnglePotential> {

        public:
            double value(vector<Coordinate> &coordinates);
    };


}


#endif //MOLFUNC_EXT_ANGLES_H
