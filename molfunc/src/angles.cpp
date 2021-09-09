#include "cmath"
#include "vector"
#include "angles.h"


using namespace std;


namespace molfunc{

    AnglePotential::AnglePotential() = default;

    AnglePotential::AnglePotential(double phi0, double k) {
        /**************************************************
         * Harmonic potential for an angle
         *
         *      V = k/2 (cos(phi) - cos(phi0))^2
         *
         * where phi is the current value of the angle and
         * the difference in the cosines is used to save a
         * arccos operation.
         *
         * Arguments:
         *      phi0: Optimal angle in degrees
         *
         *      k: Force constant
         **************************************************/

        double deg_to_rad = 0.0174533;
        this->cos_phi0 = cos(deg_to_rad * phi0);

        this->half_k = k / 2.0;
    }

    double AnglePotential::value(vector<Coordinate> &coordinates) {
        /**************************************************
         * Evaluate the value of the potential
         *
         *      V = k/2 (cos(phi) - cos(phi0))^2
         *
         * where the angle is between 0-1-2
         **************************************************/
        if (atom_idxs[0] > coordinates.size()
            || atom_idxs[0] == atom_idxs[1]
            || atom_idxs[1] > coordinates.size()
            || atom_idxs[1] == atom_idxs[2]
            || atom_idxs[2] > coordinates.size()
            || atom_idxs[0] == atom_idxs[2]){

            throw runtime_error("Invalid indexing for these coordinates");
        }

        auto vec1 = coordinates[atom_idxs[0]] - coordinates[atom_idxs[1]];
        vec1.normalise();

        auto vec2 = coordinates[atom_idxs[2]] - coordinates[atom_idxs[1]];
        vec2.normalise();

        double cos_phi = vec1.dot(vec2);
        double delta_phi = cos_phi - cos_phi0;

        return half_k * delta_phi * delta_phi;
    }


    AnglePotentials::AnglePotentials() = default;

    double AnglePotentials::value(vector<Coordinate> &coordinates) {

        double total = 0.0;

        for (auto &potential : potentials){
            total += potential.value(coordinates);
        }

        return total;
    }
}

