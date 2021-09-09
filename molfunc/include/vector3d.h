#ifndef MOLFUNC_EXT_VECTOR3D_H
#define MOLFUNC_EXT_VECTOR3D_H
#include "array"
#include "iostream"
#include "string"

using namespace std;

namespace molfunc{

    class Vector3D: public array<double, 3>{

        public:
            Vector3D operator/ (double x);
            Vector3D operator* (double x);
            Vector3D operator+ (double x);

            double length();

            void normalise();

            double dot(Vector3D &other);

    };

}

#endif //MOLFUNC_EXT_VECTOR3D_H
