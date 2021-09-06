#include "rotation.h"
#include "cmath"

using namespace std;


namespace molfunc{

    void RotationMatrix::update(double w1, double w2, double w3) {
        /*******************************************************
         *  Calculate the rotation matrix. Notation from
         *  [G. Terzakis, M. Lourakis, D. Ait-Boudaoud, J Math Imaging Vis
         *  (2018) 60:422–442 (https://doi.org/10.1007/s10851-017-0765-x)]
         *
         *     R = I_3 + sin(θ)/θ × [ω]_x + (1 - cos(θ))/θ^2 × [ω]_x^2
         *
         *  where θ = |ω|, I_3 is the 3x3 identity matrix and [ω]_x^2 is
         *  the matrix (rather than element wise) square. ω is the vector
         *  of (w1, w2, w3)
         *
         *     [ω]_x = [[0, -w3, w2],
         *              [w3, 0, -w1],
         *              [-w2, w1, 0]])
         */
        double theta = sqrt(w1*w1 + w2*w2 + w3*w3);

        double a = sin(theta) / theta;                    // sin(θ)/θ
        double b = (1.0 - cos(theta)) / (theta * theta);  // (1 - cos(θ))/θ^2

        // First row of the matrix
        this->data()[0][0] = 1 - b * (w3 * w3 + w2 * w2);
        this->data()[0][1] = - a * w3 + b * w2 * w1;
        this->data()[0][2] = a * w2 + b * w3 * w1;

        // Second row of the matrix
        this->data()[1][0] = a * w3 + b * w1 * w2;
        this->data()[1][1] = 1 - b * (w3 * w3 + w1 * w1);
        this->data()[1][2] = - a * w1 + b * w2 * w3;

        // Third row of the matrix
        this->data()[2][0] = -a * w2 + b * w1 * w3;
        this->data()[2][1] = a * w1 + b * w2 * w3;
        this->data()[2][2] = 1 - b * (w1 * w1 + w2 * w2);
    }

    void RotationMatrix::update(GridPoint &arr) {
        update(arr[0], arr[1], arr[2]);
    }
}
