import numpy as np
import scipy

from voila.constants import MINVAL
from voila.vlsv import get_expected_dpsi, matrix_area, get_expected_psi


def unpack_means(value):
    if np.size(value, 0) == 1:
        value = np.append(value, np.array(1 - value[0]))
    return value.tolist()


def unpack_bins(value):
    if np.size(value, 0) == 1:
        value = np.append(value, [np.flip(value[-1], 0)], axis=0)
    return value.tolist()


def generate_excl_incl(means):
    l = []
    for mean in means:
        if mean < 0:
            l.append((-mean, 0))
        else:
            l.append((0, mean))
    return l


def generate_means(bins):
    m = []
    for b in bins:
        m.append(get_expected_dpsi(b))
    return m


def generate_high_probability_non_changing(ir, prior, non_changing_threshold, bins):
    x = []
    prior = prior[1 if ir else 0]

    for bin in bins:
        bin = np.array(bin)
        bin += MINVAL
        bin /= bin.sum()
        A = np.log(bin) - prior
        R = np.exp(A - scipy.special.logsumexp(A))
        x.append(matrix_area(R, non_changing_threshold, non_changing=True))

    return x


def generate_variances(bins):
    v = []
    for b in bins:
        epsi = get_expected_psi(b)
        # b used to be a numpy array and now it's a list...
        step_bins = 1.0 / len(b)
        projection_prod = b * np.arange(step_bins / 2, 1, step_bins) ** 2
        v.append(np.sum(projection_prod) - epsi ** 2)
    return v
