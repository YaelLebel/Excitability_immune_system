import numpy as np

def calculate_pareto_front(points):
    """
    Calculate the Pareto front of a set of 2D points (minimization).
    Refactored from utils_pp.py.

    Parameters:
        points (array-like): A collection of points as a list of tuples or a 2D array, where each point is (x, y).

    Returns:
        pareto_indices (list): A list of indices of the points on the Pareto front in the original list.
    """
    points = np.array(points)
    # Sort by first objective
    sorted_indices = np.lexsort((points[:, 1], points[:, 0]))
    sorted_points = points[sorted_indices]

    pareto_indices = []
    current_min_y = float('inf')

    # Iterate through the sorted points
    for i, point in enumerate(sorted_points):
        x, y = point
        if y < current_min_y:
            pareto_indices.append(sorted_indices[i])
            current_min_y = y

    return pareto_indices

def calculate_pareto_front_max(points):
    """
    Calculate the Pareto front of a set of 2D points for MAXIMIZATION.
    Refactored from utils_pp.py.

    Parameters:
        points (array-like): A collection of points as a list of tuples or a 2D array, where each point is (x, y).

    Returns:
        pareto_indices (list): A list of indices of the points on the Pareto front in the original list.
    """
    points = np.array(points)
    # Negate objectives to use lexicographical sort for maximization logic essentially
    # Or just sort descending
    # lexsort sorts by last key, then second to last...
    # We want max x, max y.
    # Sort descending by x.
    
    # Sort by -y, then -x to iterate?
    # utils_pp logic: sorted_indices = np.lexsort((-points[:, 1], -points[:, 0]))
    
    sorted_indices = np.lexsort((-points[:, 1], -points[:, 0]))
    sorted_points = points[sorted_indices]

    pareto_indices = []
    current_max_y = -float('inf')

    for i, point in enumerate(sorted_points):
        x, y = point
        # If we sorted by X descending, then for a point to be on pareto front (max-max),
        # its Y must be greater than any max Y seen so far?
        # Standard sweep line algorithm.
        if y > current_max_y:
            pareto_indices.append(sorted_indices[i])
            current_max_y = y

    return pareto_indices
