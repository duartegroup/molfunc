import numpy as np


def get_rotated_coords(rotation, coords):
    """Rotate coordinates by theta radians in an axis

    :param rotation: (np.ndarray) shape: (4,). [x, y, z, theta]
    :param coords: (np.ndarray) shape: (n_atoms, 3)
    :return: rotated coordinates
    """
    rot_mat = rotation_matrix(axis=rotation[:3], theta=rotation[3])

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


def get_points_on_sphere(n_points, r=1):
    """
    Find n evenly spaced points on a sphere using the "How to generate
    equidistributed points on the surface of a sphere" by Markus Deserno, 2004.

    Arguments:
        n_points (int): number of points to generate
        r (float): radius of the sphere

    Returns:
        (list(np.ndarray))
    """
    points = []

    a = 4.0 * np.pi * r**2 / n_points
    d = np.sqrt(a)
    m_theta = int(np.round(np.pi / d))
    d_theta = np.pi / m_theta
    d_phi = a / d_theta

    for m in range(m_theta):
        theta = np.pi * (m + 0.5) / m_theta
        m_phi = int(np.round(2.0 * np.pi * np.sin(theta)/d_phi))

        for n in range(m_phi):
            phi = 2.0 * np.pi * n / m_phi
            point = np.array([r * np.sin(theta) * np.cos(phi),
                              r * np.sin(theta) * np.sin(phi),
                              r * np.cos(theta)])

            points.append(point)

    # This algorithm doesn't always generate enough points (at least not quite)
    # so add a unit vector on the x axis to make up the difference
    if len(points) < n_points:
        points += (n_points - len(points)) * [np.array([1.0, 0.0, 0.0])]

    return points

