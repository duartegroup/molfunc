#include "vector3d.h"


namespace molfunc{

    Vector3D Vector3D::operator*(double x) {

        return {this->data()[0] * x,
                this->data()[1] * x,
                this->data()[2] * x};
    }

    Vector3D Vector3D::operator/(double x) {

        return {this->data()[0] / x,
                this->data()[1] / x,
                this->data()[2] / x};
    }

    Vector3D Vector3D::operator+(double x) {

        return {this->data()[0] + x,
                this->data()[1] + x,
                this->data()[2] + x};
    }
}

