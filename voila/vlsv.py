import numpy as np


def get_expected_value(bins, left=0, right=1):
    """
    Get numpy array of expected values.
    :param bins: bins from majiq
    :param left: left constraint
    :param right: right constraint
    :return: numpy array
    """

    step = (right - left) / len(bins)
    return np.arange(left + step / 2, right, step).dot(bins)


def get_expected_dpsi(bins):
    """
    Get numpy array of expected delta psi values.
    :param bins: bins from majiq
    :return: numpy array
    """

    return get_expected_value(bins, left=-1)


def get_expected_psi(bins):
    """
    Get numpy array of expect psi values.
    :param bins: bins from majiq
    :return: numpy array
    """

    return get_expected_value(bins)


def collapse_matrix(matrix):
    """
    Collapse the diagonals probabilities in 1-D and return them.
    :param matrix: numpy matrix
    :return: collapsed numpy array
    """

    collapse = []
    matrix_corner = matrix.shape[0]
    for i in range(-matrix_corner + 1, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)


def matrix_area(matrix, threshold, non_changing=False):
    """
    Returns the probability of an event to be above a certain threshold.  If non_changing is set, then return area for
    non changing events.
    :param non_changing: boolean
    :param threshold: psi threshold
    :param matrix: numpy matrix
    :return: probability
    """

    collapse = matrix
    collapse = np.concatenate(([0], collapse))
    collapse = np.cumsum(collapse)
    collapse /= collapse[-1]
    xbins = np.linspace(-1, 1, num=collapse.size)
    abs_threshold = abs(threshold)
    left, right = np.interp([-abs_threshold, abs_threshold], xbins, collapse, left=0, right=1)
    if non_changing:
        area = right - left
    else:
        area = left + (1 - right)
    return np.clip(area, 0, 1)


def is_lsv_changing(means, threshold):
    """
    Return true if lsv is changing based on threshold.
    :param threshold: lsv threshold value
    :return: bool
    """
    means = np.array(tuple(means))
    means_gt_zero = means[means > 0]
    means_sum = means_gt_zero.sum()
    max_value = max(means_sum, abs(means_sum))
    return max_value >= threshold
