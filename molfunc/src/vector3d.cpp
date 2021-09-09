#include "cmath"
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

    void Vector3D::normalise(){
        /**********************
         * Normalise the vector
         **********************/

        double l = length();

        this->data()[0] /= l;
        this->data()[1] /= l;
        this->data()[2] /= l;
    }

    double Vector3D::length() {
        /**********************
         * L2 norm
         **********************/

        double length_sq = this->data()[0] * this->data()[0]
                           + this->data()[1] * this->data()[1]
                           + this->data()[2] * this->data()[2];

        return sqrt(length_sq);
    }

    double Vector3D::dot(Vector3D &other) {
        /****************************
         * Dot product of two vectors
         ****************************/

        return this->data()[0] * other[0]
               + this->data()[1] * other[1]
               + this->data()[2] * other[2];
    }
}

