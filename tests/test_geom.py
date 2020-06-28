from molfunc_ext import get_rotated_coordinates
import numpy as np


def test_molecule_rotation():

    coords = np.array([[1.500e-03,  2.000e-04,  3.400e-03],
                       [-5.948e-01, -2.155e-01,  9.059e-01],
                       [4.716e-01,  9.986e-01,  2.830e-02],
                       [-6.780e-01, -3.370e-02, -8.588e-01],
                       [7.999e-01, -7.496e-01, -7.890e-02]])

    rot_coords = [[6.42159262e-04,  1.37026694e-03,  3.40000000e-03],
                  [-1.40034814e-01, -6.16942089e-01,  9.05900000e-01],
                  [-5.85486358e-01,  9.36383599e-01,  2.83000000e-02],
                  [-3.37967391e-01, -5.88725515e-01, -8.58800000e-01],
                  [1.06295446e+00,  2.68082032e-01, -7.89000000e-02]]

    new_coords = get_rotated_coordinates(py_omega=[0.0, 0.0, 1.0],
                                         py_coords=coords.flatten())

    # Check that the new coordinates are close to the expected (rot_coords)
    for i, coord in enumerate(new_coords):
        assert np.linalg.norm(coord - rot_coords[i]) < 1E-6

