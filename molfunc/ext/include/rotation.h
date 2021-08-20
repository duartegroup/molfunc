#ifndef MOLFUNC_EXT_ROTATION_H
#define MOLFUNC_EXT_ROTATION_H
#include "array"


using namespace std;

namespace molfunc {
    class RotationMatrix : public array<array<double, 3>, 3> {

    public:

        void update(double w1, double w2, double w3);
        void update(array<double, 3> &arr);

    };

}

#endif //MOLFUNC_EXT_ROTATION_H
