from molfunc.geom import get_rotated_coords
from scipy.spatial import distance_matrix
from scipy.optimize import minimize
from molfunc.geom import get_points_on_sphere
import numpy as np


def repulsion(rotation, coords, alt_coords):
    """
    Calculate the repulsion energy between two sets of coordinates with

        V = Σpairs 1/rij^4

    where the coordinates are rotated in an axis x[:3] theta radians (x[3])

    :param rotation: (np.ndarray) shape: (4,)
    :param coords: (np.ndarray) shape: (n_atoms, 3)
    :param alt_coords: (np.ndarray) shape: (n_other_atoms, 3)
    """

    coords = get_rotated_coords(rotation, coords)
    dist_mat = distance_matrix(coords, alt_coords)
    energy = np.sum(np.power(dist_mat, -4))

    return energy


def get_lowest_energy_rotation(coords, other_coords, n=100):
    """
    Get the best rotation array given a set of coords, that can be rotated to
    minimise the energy with respect to another set of coordinates

    :param coords: (np.ndarray) shape (n, 3)
    :param other_coords: (np.ndarray) shape (m, 3)
    :param n: (int) Number of trials to (hopefully) find the global minimum

    :return: (np.ndarray) shape: (4,) Rotation array of an axis (length 3)
             and an angle (rad) to be used in get_rotated_coords()
    """
    # TODO parallelise
    min_energy = best_rotation = None

    # Get a set of equally distributed points on a sphere to use as axes
    axes = get_points_on_sphere(n_points=n)
    rotations = [np.append(axes[i], np.random.uniform(0.0, 2*np.pi)) for i in range(n)]

    for rotation in rotations:

        # Minimise repulsion between coords and the other_coords
        # TODO analytic derivative
        res = minimize(repulsion, x0=rotation, args=(coords, other_coords),
                       method='L-BFGS-B', tol=1.0)

        # Minimise the energy to a tighter tolerance if it's fairly low energy
        if res.fun < 5:
            res = minimize(repulsion, x0=res.x, args=(coords, other_coords),
                           method='L-BFGS-B', tol=0.05)

        if min_energy is None or res.fun < min_energy:
            min_energy = res.fun
            best_rotation = res.x

    return best_rotation
