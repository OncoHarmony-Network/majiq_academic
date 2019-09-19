import numpy as np
from scipy import special

from voila.constants import MINVAL
from voila.vlsv import get_expected_dpsi, matrix_area, get_expected_psi
from voila.voila_log import voila_log


def unpack_means(value):
    """
    In the case where lsv is binary, we need to generate the other junction data.
    :param value: means
    :return: list
    """
    if np.size(value, 0) == 1:
        value = np.append(value, np.array(1 - value[0]))
    return value.tolist()


def unpack_bins(value):
    """
    In the case where lsv is binary, we need to generate the other junction data.
    :param value: bins
    :return: list
    """
    if np.size(value, 0) == 1:
        value = np.append(value, [np.flip(value[-1], 0)], axis=0)
    return value.tolist()


def generate_excl_incl(means):
    """
    Create exclusion and inclusion values for plots.
    :param means: lsv means data
    :return: list
    """
    l = []
    for mean in means:
        if mean < 0:
            l.append((-mean, 0))
        else:
            l.append((0, mean))
    return l


def generate_means(bins):
    """
    Create means where not available in voila file.
    :param bins: bins data
    :return: list
    """
    m = []
    for b in bins:
        m.append(get_expected_dpsi(b))
    return m


def generate_prior_removed_expected_dpsi(ir, prior, bins=None, pr_removed_bins=None):
    """
    Similar input / output use case to the function below, but gets expected dpsi instead of matrix area (prob)
    """
    x = []

    for pr_removed_bin in generate_bins_prior_removed(ir, prior, bins) if bins is not None else pr_removed_bins:
        x.append(abs(get_expected_dpsi(pr_removed_bin)))

    return x

def generate_high_probability_non_changing(ir, prior, non_changing_threshold, bins=None, pr_removed_bins=None):
    """
    Calculate the probability of non changing lsv junctions.
    :param ir: Does this lsv have intron retention.
    :param prior: prior matrix from voila file.
    :param non_changing_threshold: non-changing threshold set by user.
    :param bins: bins data from voila file.
    :return: list
    """
    x = []

    for pr_removed_bin in generate_bins_prior_removed(ir, prior, bins) if bins is not None else pr_removed_bins:
        x.append(matrix_area(pr_removed_bin, non_changing_threshold, non_changing=True))

    return x

def generate_bins_prior_removed(ir, prior, bins):
    new_bins = []
    prior = prior[1 if ir else 0]

    for bin in bins:
        bin = np.array(bin)
        bin += MINVAL
        bin /= bin.sum()
        A = np.log(bin) - prior
        R = np.exp(A - special.logsumexp(A))
        new_bins.append(R)

    return new_bins



def generate_variances(bins):
    """
    Calulate variances for lsv junctions.
    :param bins: bins data from voila file.
    :return: list
    """
    v = []
    for b in bins:
        epsi = get_expected_psi(b)
        # b used to be a numpy array and now it's a list...
        step_bins = 1.0 / len(b)
        projection_prod = b * np.arange(step_bins / 2, 1, step_bins) ** 2
        var = np.sum(projection_prod) - (epsi ** 2)
        if var < 0:
            v.append(0.0)
            if abs(var) < step_bins:
                # if var is < 0 but also of less magnitude then the bin size, we assume a rounding error and
                # pass along a zero variance
                pass
            else:
                # otherwise, this should not happen, so we send a warning
                log = voila_log()
                log.warning("Negative variance found when calculating over %s" % str(b))
        else:
            v.append(var)


    return v
