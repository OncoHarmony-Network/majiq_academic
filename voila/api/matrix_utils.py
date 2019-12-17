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


def generate_high_probability_non_changing(ir, prior, non_changing_threshold, bins):
    """
    Calculate the probability of non changing lsv junctions.

    Done with respect to posterior with prior factored out to focus on 'truly
    non-changing events' that are not non-changing only because of the prior.
    :param ir: Does this lsv have intron retention.
    :param prior: prior matrix from voila file (absolute probabilities).
    :param non_changing_threshold: non-changing threshold set by user.
    :param bins: bins data from voila file (absolute probabilities).
    :return: list
    """
    x = []
    prior = prior[1 if ir else 0]

    for bin in bins:
        bin = np.array(bin)
        bin += MINVAL
        bin /= bin.sum()
        # factor out prior from posterior (in logspace to reduce underflow)
        A = np.log(bin) - np.log(prior)
        # subtract logsumexp to renormalize to 1, then transform back
        R = np.exp(A - special.logsumexp(A))
        x.append(matrix_area(R, non_changing_threshold, non_changing=True))

    return x

def generate_standard_deviations(bins):
    """
    Calulate standard deviations for lsv junctions.
    :param bins: bins data from voila file.
    :return: list
    """
    v = np.sqrt(generate_variances(bins))
    return v

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
