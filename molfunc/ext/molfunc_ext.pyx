# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
from cpython.array cimport array, clone
from libc.math cimport sqrt, pow, cos, sin
import numpy as np


cdef get_rotation_matrix(float w1, float w2, float w3):
    """
    Calculate the rotation matrix. Notation from [G. Terzakis, M. Lourakis,
    D. Ait-Boudaoud, J Math Imaging Vis (2018) 60:422–442
    (https://doi.org/10.1007/s10851-017-0765-x)]

    R = I_3 + sin(θ)/θ × [ω]_x + (1 - cos(θ))/θ^2 × [ω]_x^2

    where θ = |ω|, I_3 is the 3x3 identity matrix and [ω]_x^2 is the matrix
    (rather than element wise) square. ω is the vector of (w1, w2, w3)

    [ω]_x = [[0, -w3, w2],
             [w3, 0, -w1],
             [-w2, w1, 0]])

    :param w1: (float)
    :param w2: (float)
    :param w3: (float)
    :return: (array) Rotation matrix (R). shape = (3, 3)
    """
    cdef double R[3][3]

    # θ = |ω|
    cdef double theta = sqrt(w1**2 + w2**2 + w3**2)

    # sin and cos constants
    cdef double sin_factor = sin(theta) / theta               # sin(θ)/θ
    cdef double cos_factor = (1.0 - cos(theta)) / theta**2    # (1 - cos(θ))/θ^2

    # First row of the matrix
    R[0][0] = 1 - cos_factor * (w3**2 + w2**2)
    R[0][1] = - sin_factor * w3 + cos_factor * w2 * w1
    R[0][2] = sin_factor * w2 + cos_factor * w3 * w1

    # Second row of the matrix
    R[1][0] = sin_factor * w3 + cos_factor * w1 * w2
    R[1][1] = 1 - cos_factor * (w3**2 + w1**2)
    R[1][2] = - sin_factor * w1 + cos_factor * w2 * w3

    # Thrid row of the matrix
    R[2][0] = -sin_factor * w2 + cos_factor * w1 * w3
    R[2][1] = sin_factor * w1 + cos_factor * w2 * w3
    R[2][2] = 1 - cos_factor * (w1**2 + w2**2)

    return R


cdef rotate_coordinates(array coords, int n_atoms, array omega):
    """
    Rotate coordinates using a rotation matrix given a vector omega

    :param coords: (array) length = 3 * n_atoms
    :param n_atoms: (int)
    :param omega: (array) length = 3
    """
    cdef double R[3][3]
    R = get_rotation_matrix(omega[0], omega[1], omega[2])

    cdef int i
    cdef double x, y, z

    # Apply the roation matrix to the flat array of coordinates
    for i in range(n_atoms):

        # Tempoary values used to compute the matrix multiplication
        x = coords.data.as_doubles[3*i]
        y = coords.data.as_doubles[3*i+1]
        z = coords.data.as_doubles[3*i+2]

        # x_i
        coords.data.as_doubles[3*i] = (R[0][0] * x + R[0][1] * y + R[0][2] * z)

        # y_i
        coords.data.as_doubles[3*i+1] = (R[1][0] * x + R[1][1] * y + R[1][2] * z)

        # z_i
        coords.data.as_doubles[3*i+2] = (R[2][0] * x + R[2][1] * y + R[2][2] * z)

    return None


cdef get_energy(array coords_i, int n_atoms_i, array coords_j, int n_atoms_j):
    """
    Given two sets of coordinates calculate the replsive energy between them as

    energy = Σpairs 1/rij^8

    :param coords_i: (array) length = 3 * n_atoms_i e.g.
                     [x_1, y_1, z_1, x_2, y_2, z_2, ....]
    :param coords_j: (array) length = 3 * n_atoms_j
    :return: (float) energy
    """
    cdef int i, j
    cdef float energy = 0.0

    cdef double delta_x
    cdef double delta_y
    cdef double delta_z
    cdef double d

    for i in range(n_atoms_i):
        for j in range(n_atoms_j):

            # Distance between coordinate i and coordinate j
            delta_x = coords_i.data.as_doubles[3*i] - coords_j.data.as_doubles[3*j]
            delta_y = coords_i.data.as_doubles[3*i+1] - coords_j.data.as_doubles[3*j+1]
            delta_z = coords_i.data.as_doubles[3*i+2] - coords_j.data.as_doubles[3*j+2]
            d_sq = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z

            energy += pow(d_sq, -4)

    return energy


cpdef get_array(py_array):
    """
    Convert a Python list or numpy ndarray into a C array. Size (n) is
    unknown at compile time so needs templates

    :param py_array: (list, np.ndarray)
    """

    cdef array template = array('d')

    cdef int n = int(len(py_array))
    cdef c_array = clone(template, n, False)

    # Populate the array
    cdef int i
    for i in range(n):
        c_array[i] = py_array[i]

    return c_array


cpdef get_rotated_coordinates(py_coords, py_omega):
    """From an ω vector rotate a set of flat coordinates -> n_atoms x 3"""

    # Set the C arrays
    cdef array coords = get_array(py_coords)
    cdef array omega = get_array(py_omega)

    # Do the rotation
    cdef int n_atoms = int(len(py_coords) / 3)
    rotate_coordinates(coords, n_atoms, omega)

    # Return the n_atoms x 3 numpy array
    return np.array(coords).reshape(n_atoms, 3)


cpdef get_minimised_coords(py_coords, py_other_coords):
    """
    Get the best rotation array given a set of coords, that can be rotated to
    minimise the energy with respect to another set of coordinates.

    :param coords: (np.ndarray) shape (n, 3)
    :param other_coords: (np.ndarray) shape (m, 3)
    """
    cdef int n_atoms_i = len(py_coords)
    cdef int n_atoms_j = len(py_other_coords)

    # Faster to work with flat arrays
    py_coords, py_other_coords = py_coords.flatten(), py_other_coords.flatten()

    # Convert to C arrays
    cdef array coords_i = get_array(py_coords)
    cdef array coords_j = get_array(py_other_coords)

    # Have a copy of the coordinates to reset following rotation and caluclating
    # the energy
    cdef array coords_i_copy = get_array(py_coords)

    # Grid values for the rotation-matrix defining rotation
    cdef double[12] w1s = np.linspace(-100, 100, 12)
    cdef double[12] w2s = np.linspace(-100, 100, 12)
    cdef double[12] w3s = np.linspace(-100, 100, 12)

    # Iterators
    cdef int i, j, k

    # Minimum energy and the best ω vector (that generates the lowest energy)
    cdef double min_energy = 9999999.9
    cdef double energy = 0.0

    # current and optimal ω vectors
    cdef array omega = get_array(np.ones(3))
    cdef array best_omega = get_array(np.ones(3))

    # Enumerate all possible rotation values
    for i in range(12):
        for j in range(12):
            for k in range(12):
                # Set the ω vector
                omega[0] = w1s[i]
                omega[1] = w2s[j]
                omega[2] = w3s[k]

                # Rotate atoms in set i with omega in place
                rotate_coordinates(coords_i, n_atoms_i, omega)

                # Caluclate the repulsive energy between the coordinates
                energy = get_energy(coords_i, n_atoms_i, coords_j, n_atoms_j)

                # If the energy is lower then reset best_omega
                if energy < min_energy:
                    min_energy = energy
                    best_omega[0] = w1s[i]
                    best_omega[1] = w2s[j]
                    best_omega[2] = w3s[k]

                # Reset the values of coordinates i
                for n in range(3*n_atoms_i):
                    coords_i[n] = coords_i_copy[n]

    return get_rotated_coordinates(py_coords, best_omega)
