from molfunc.geom import rotation_matrix
import numpy as np


def test_rotation():

    point = np.array([1.0, 0.0, 0.0])

    # z axis, length 2 to check that the axis is normalised in the function
    axis = 2 * np.array([0.0, 0.0, 1.0])

    rot_matrix = rotation_matrix(axis, theta=np.pi)

    new_point = np.matmul(rot_matrix, point)
    # Rotation of 180 degrees in the z axis a point x should -> -x

    assert np.linalg.norm(new_point - np.array([-1.0, 0.0, 0.0])) < 1E-6
