import numpy as np


def get_rotated_coords(x, coords):
    """Rotate coordinates by theta radians in an axis

    :param x: (np.ndarray) shape: (4,). [x, y, z, theta]
    :param coords: (np.ndarray) shape: (n_atoms, 3)
    :return: rotated coordinates
    """
    rot_mat = rotation_matrix(axis=x[:3], theta=x[3])

    def rotate(coord):
        return np.matmul(rot_mat, coord)

    # Map the matrix multiplication over all atomic coordinates
    return np.array(list(map(rotate, coords)))


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    :param axis: (np.ndarray) Unit vector in 3D to rotate around
    :param theta: (float) Angle in radians
    """
    axis = np.asarray(axis)
    axis = axis/np.linalg.norm(axis)
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
