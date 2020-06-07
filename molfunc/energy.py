from molfunc.geom import get_rotated_coords
from scipy.spatial import distance_matrix
import numpy as np


def repulsion(x, coords, alt_coords):
    """
    Calculate the repulsion energy between two sets of coordinates with

        V = Σpairs 1/rij^4

    where the coordinates are rotated in an axis x[:3] theta radians (x[3])

    :param x: (np.ndarray) shape: (4,)
    :param coords: (np.ndarray) shape: (n_atoms, 3)
    :param alt_coords: (np.ndarray) shape: (n_other_atoms, 3)
    """

    coords = get_rotated_coords(x, coords)
    dist_mat = distance_matrix(coords, alt_coords)
    energy = np.sum(np.power(dist_mat, -4))

    return energy
