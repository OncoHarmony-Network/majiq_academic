import numpy as np


# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


def random_merge(mat_a,
                 mat_b,
                 preserve_order=False,
                 seed=False):
    """
    Given two numpy 2D arrays with same shape, return matrix with same shape
        that has randomly selected rows from first two matrices such that
        all 0->n_rows are represented in return matrix.

        Credit to Scott Norton for writing this

    :param mat_a: np.array
    :param mat_b: np.array
    :param preserve_order: Boolean, final result must have rows in same order as
        original matrices?
    :param seed: if set to an integer, makes random behavior reproducible lol
    :return: np.array
    """
    assert mat_a.shape == mat_b.shape
    if isinstance(seed, int):
        np.random.seed(seed)
    elif not isinstance(seed, bool):
        raise ValueError("seed needs to be a bool or int.")
    elif seed:  # else True or False, if True, set seed to 0
        np.random.seed(0)
    # get array of indices 0->len(array)
    mat_c = np.zeros((mat_a.shape[0], mat_a.shape[1]))
    # Permute the row indices.  This ensures that in the forthcoming loop,
    # each row index is used once and only once.
    indices = np.random.permutation(mat_a.shape[0])
    # Enter the main loop.
    if preserve_order:
        # Initialize the target matrix
        mat_c = np.zeros((mat_a.shape[0], mat_a.shape[1]))
        # Enter the main loop.
        for i in range(mat_a.shape[0]):  # row contains the element of indices,
                                           # i contains its position
            mat_c[i] = (mat_a, mat_b)[np.random.rand() < 0.5][i].copy()
            # We use the copy method here so that we can modify mat_a and mat_b later.
        return mat_c
    for i, row in enumerate(indices):  # row contains the element of indices,
                                       # i contains its position
        mat_c[i] = (mat_a, mat_b)[np.random.rand() < 0.5][row].copy()
        # We use the copy method here so that we can modify mat_a and mat_b later.
    return mat_c
