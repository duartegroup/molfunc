#include "molecules.h"
#include "utils.h"
#include "iostream"
#include "catch2/catch.hpp"


using namespace std;
using namespace molfunc;


TEST_CASE("Test H atom rotation"){

    vector<Atom3D> atoms = {Atom3D("H", 1.0, 0.0, 0.0)};
    auto mol = Molecule(atoms);

    RotationMatrix rot_mat = {{
                     {{-1.0, 0.0, 0.0},
                         {0.0, -1.0, 0.0},
                         {{0.0, 0.0, 1.0}}}
                              }};

    // Rotating by π/2 around the z axis should leave the atom at (-1, 0, 0)
    mol.rotate(rot_mat);
    REQUIRE(utils::is_close(mol.coordinates[0], {-1, 0, 0}));

    // translate back to its original position
    mol.translate({2.0, 0.0, 0.0});

    // Rotating about the y-axis should give the same result
    rot_mat = {{
               {{-1.0, 0.0, 0.0},
           {0.0, 1.0, 0.0},
           {{0.0, 0.0, -1.0}}}
                }};
    mol.rotate(rot_mat);
    REQUIRE(utils::is_close(mol.coordinates[0], {-1, 0, 0}));

}


TEST_CASE("Test rotation preserves norm for H2"){

    vector<Atom3D> atoms = {Atom3D("H", 0.0, 0.0, 0.0),
                            Atom3D("H", 1.0, 0.0, 0.0)};
    auto mol = Molecule(atoms);

    RotationMatrix rot_mat = {};  // Zero construction of 3x3 rotation matrix

    for (int i=0; i<6; i++){
        for (int j=0; j<6; j++){
            for (int k=0; k<6; k++){
                rot_mat.update(-1 + 2*static_cast<double>(i)/5,  // -1-->1
                               -1 + 2*static_cast<double>(j)/5,
                               -1 + 2*static_cast<double>(k)/5);

                mol.rotate(rot_mat);

                // Rotation should preserve the distance
                REQUIRE(utils::is_close(mol.distance(0, 1), 1.0));
            }// k
        }// j
    }// i
}


TEST_CASE("Test equal distro of rotations"){

    vector<Atom3D> atoms = {Atom3D("H", 0.0, 0.0, 0.0),
                            Atom3D("H", 1.0, 0.0, 0.0)};
    auto mol = Molecule(atoms);

    RotationMatrix rot_mat = {};     // Zero construction of 3x3 rotation matrix

    array<int, 19> histogram = {};   // Histogram of cosθ in bins of
                                     // {(-1, -0.9), (-0.9, -0.8) ...}

    int num = 10;
    auto num_d = double(num);
    double min_w = -100.0;
    double max_w = 100.0;

    for (int i=0; i<num; i++){
        for (int j=0; j<num; j++){
            for (int k=0; k<num; k++){

                rot_mat.update(min_w + (max_w - min_w)*static_cast<double>(i)/(num_d-1),  // -1-->1
                               min_w + (max_w - min_w)*static_cast<double>(j)/(num_d-1),
                               min_w + (max_w - min_w)*static_cast<double>(k)/(num_d-1));

                mol.rotate(rot_mat);

                //                   cos(angle to x axis)
                int bin_idx = int(((mol.coordinates[1].x() + 1.0)/2.0) * 20);
                histogram[bin_idx] += 1;

                mol.coordinates[1] = {1.0, 0.0, 0.0};
            }// k
        }// j
    }// i

    double total = pow(num_d, 3);
    double ideal_fraction = 19 / total;

    for (auto frequency: histogram){
        REQUIRE(utils::is_close(frequency/total, ideal_fraction, 0.2));

    }
}
