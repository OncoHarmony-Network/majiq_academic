import random
from numpy.ma import masked_less
import numpy as np
from scipy.stats import nbinom, poisson
# import majiq.src.plotting as mplot

def get_ecdf(pvalues):
    # print sorted(pvalues)
    nbins = max(min(10, len(pvalues)), len(pvalues) / 10)
    hist, bin_edges = np.histogram(pvalues, range=[0, 1], bins=nbins, density=True)
    return np.cumsum(hist) / len(bin_edges)


def score_ecdf(ecdf):
    """
    Give a score to a ecdf calculating the deviation from the 45 degree line
    """
    return sum(abs(np.linspace(0, 1, num=len(ecdf)) - ecdf))


def __calc_nbin_p(r, mu):
    p = r / (r + mu)
    return p


def sample_over_nb(one_over_r, mu, num_samples):
    if one_over_r > 0:
        r = 1 / one_over_r
        p = __calc_nbin_p(r, mu)
        sampl = nbinom.rvs(r, p, size=num_samples)
    else:
        sampl = poisson.rvs(mu, size=num_samples)
    return sampl


def get_negbinom_pval(one_over_r, mu, x):
    if one_over_r > 0:
        r = 1 / one_over_r
        p = __calc_nbin_p(r, mu)
        nbcdf = nbinom.cdf(x, r, p)  # + nbinom.pmf(x, r, p)
    else:
        nbcdf = poisson.cdf(x, mu)
    return 1 - nbcdf


def calc_pvalues(junctions, one_over_r, indices_list=None):
    pvalues = []
    for i, junc in enumerate(junctions):

        # get mu and jpos
        if indices_list is None:
            junc = junc[junc.nonzero()]
            jpos = random.choice(junc)
        else:
            jpos = junc[indices_list[i]]

        ljunc = len(junc.nonzero()[0])
        mu = float(junc.sum() - jpos) / float(ljunc - 1)
        pval = get_negbinom_pval(one_over_r, mu, jpos)
        # pval = get_ztnbin_pval(one_over_r, mu, jpos)
        pvalues.append(pval)

    return pvalues


def adjust_fit(starting_a, junctions, precision, previous_score, plotpath, indices=None, logger=None):
    previous_a = -1
    if logger:
        logger.debug("Starting from %s with precision %s" % (starting_a, precision))
    idx = 0

    steps = np.arange(starting_a, 0, -precision)
    steps = np.append(steps, 0)
    for corrected_a in steps:

        # since we are reducing the "a" from the fit and the problem is too much variability, we
        # expect optimization to be getting the "a" below

        pvalues = calc_pvalues(junctions, corrected_a, indices)
        ecdf = get_ecdf(pvalues)
        score = score_ecdf(ecdf)
        # mplot.plot_fitting(ecdf, plotpath, title="%s.[step %d] 1\_r %s" % (precision, idx, corrected_a),
        #                    plotname='%s.step%s' % (precision, idx))
        idx += 1
        if logger:
            logger.debug("New Score %.5f" % score)
        if previous_score < score:
            # the best fit are previous_a and previous_score
            if previous_a == -1:
                return corrected_a, score, ecdf, pvalues
            else:
                return previous_a, previous_score, previous_ecdf, previous_pvalues

        elif corrected_a == 0:
            return corrected_a, score, ecdf, pvalues

        previous_a = corrected_a
        previous_score = score
        previous_ecdf = ecdf
        previous_pvalues = pvalues
        # pvalues = []

    if logger:
        logger.warning("WARNING: Something is wrong, please contact Biociphers!")
    return corrected_a, score, ecdf, pvalues
    # this return should not be hit


def fit_nb(junctionl, outpath, plotpath, nbdisp=0.1, logger=None):
    if logger and plotpath:
        logger.debug("NBFit: Plots will be drawn in %s..." % plotpath)

    filtered = []

    for jdx, jun in enumerate(junctionl):
        if np.count_nonzero(jun) >= 5 and jun.sum() >= 10:
            filtered.append(jun)

    junctions = np.array(filtered)
    junctions = masked_less(junctions, 0.1)
    mean_junc = junctions.mean(axis=1)
    std_junc = junctions.std(axis=1)

    indices = np.zeros(shape=len(junctions), dtype=np.int)
    for i, jj in enumerate(junctions):
        jji = jj.nonzero()
        indices[i] = np.random.choice(jji[0])

    # linear regression, retrieve the a and the b plus
    one_over_r0, b = np.polyfit(mean_junc, std_junc, 1)

    pvalues = calc_pvalues(junctions, one_over_r0, indices)
    ecdf = get_ecdf(pvalues)
    # mplot.plot_fitting(ecdf, plotpath, title="NON-Corrected ECDF 1\_r %s" % one_over_r0)
    # plot_negbinomial_fit(mean_junc, std_junc, fit_function, plotpath, "Before correction")
    score = score_ecdf(ecdf)

    precision_values = [0.1, 0.01, 0.001]

    one_over_r = one_over_r0

    for i, precision in enumerate(precision_values):
        one_over_r, score, ecdf, pvalues = adjust_fit(one_over_r, junctions, precision, score, plotpath,
                                                      indices=indices, logger=logger)
        if logger:
            logger.debug("Corrected to %.5f with precision %s. Current score is %.5f" % (one_over_r, precision, score))
        if i + 1 != len(precision_values):
            #go "up" in the scale so we dont miss better solution
            one_over_r += precision - precision_values[i + 1]
            pvalues = calc_pvalues(junctions, one_over_r, indices)
            ecdf = get_ecdf(pvalues)
            score = score_ecdf(ecdf)

    if logger:
        logger.debug("Calculating the nb_r and nb_p with the new fitted function")

    return one_over_r