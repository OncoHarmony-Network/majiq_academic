import random
from numpy.ma import masked_less
import numpy as np
cimport numpy as np
from scipy.stats import nbinom, poisson
import cython

ctypedef np.float64_t DTYPE_t

# import majiq.src.plotting as mplot

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray[np.float, ndim=1] _get_ecdf(np.ndarray[DTYPE_t, ndim=1] pvalues):
    cdef int nbins
    cdef np.ndarray hist, bin_edges

    # print sorted(pvalues)
    nbins = max(min(10, len(pvalues)), len(pvalues) / 10)
    hist, bin_edges = np.histogram(pvalues, range=[0, 1], bins=nbins, density=True)
    return np.cumsum(hist) / len(bin_edges)

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef float _score_ecdf(np.ndarray[DTYPE_t, ndim=1] ecdf):
    """
    Give a score to a ecdf calculating the deviation from the 45 degree line
    """

    return sum(abs(np.linspace(0, 1, num=len(ecdf)) - ecdf))

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef inline float __calc_nbin_p(float r, float mu):
    p = r / (r + mu)
    return p

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef float get_negbinom_pval(float one_over_r, float mu, float x):
    if one_over_r > 0:
        r = 1 / one_over_r
        p = __calc_nbin_p(r, mu)
        nbcdf = nbinom.cdf(x, r, p)  # + nbinom.pmf(x, r, p)
    else:
        nbcdf = poisson.cdf(x, mu)
    return 1 - nbcdf

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray _calc_pvalues(np.ndarray[DTYPE_t, ndim=2] junctions, float one_over_r,
                              np.ndarray[np.int_t, ndim=1] indices_list):
    cdef np.ndarray[DTYPE_t, ndim=1] pvalues,
    cdef int njuncs, idx
    cdef np.ndarray[DTYPE_t, ndim=1] junc_idxs, mu, xx
    cdef np.ndarray[DTYPE_t, ndim=2] junc_fltr
    cdef np.ndarray[np.int_t, ndim=1] vals

    vals = np.count_nonzero(junctions, axis=1)
    junc_fltr = junctions[vals > 1]
    junc_idxs = np.array([xx[indices_list[idx]] for idx, xx in enumerate(junc_fltr)])
    vals = vals[vals > 1] - 1
    njuncs = vals.shape[0]

    mu = (junc_fltr.sum(axis=1) - junc_idxs) / vals
    if one_over_r > 0:
        r = 1 / one_over_r
        p = r/ (mu +r)
        pvalues = nbinom.cdf(junc_idxs, r, p)
    else:
        pvalues = poisson.cdf(junc_idxs, mu)

    return 1 - pvalues


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef list _calc_pvalues_old(np.ndarray[DTYPE_t, ndim=2] junctions, float one_over_r,
                            np.ndarray[np.int_t, ndim=1] indices_list):
    cdef list pvalues,
    cdef int njuncs, idx
    cdef np.ndarray junc_fltr, junc_idxs, vals, mu, xx

    vals = np.count_nonzero(junctions, axis=1)
    junc_fltr = junctions[vals > 1]
    junc_idxs = np.array([xx[indices_list[idx]] for idx, xx in enumerate(junc_fltr)])
    vals = vals[vals > 1] - 1
    njuncs = vals.shape[0]

    mu = (junc_fltr.sum(axis=1) - junc_idxs) / vals
    pvalues = [get_negbinom_pval(one_over_r, mu[idx], junc_idxs[idx]) for idx in range(njuncs)]

    return pvalues




@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple _adjust_fit(float starting_a, np.ndarray[DTYPE_t, ndim=2] junctions, float precision,
                       float previous_score, str plotpath, np.ndarray[np.int_t, ndim=1] index,
                       object logger=None):


    cdef int previous_a = -1
    cdef int idx = 0
    cdef np.ndarray pvalues, previous_pvalues, steps = np.arange(starting_a, 0, - precision)
    cdef np.ndarray ecdf, previous_ecdf
    cdef float score

    if logger:
        logger.debug("Starting from %s with precision %s" % (starting_a, precision))

    steps = np.append(steps, 0)
    for corrected_a in steps:

        # since we are reducing the "a" from the fit and the problem is too much variability, we
        # expect optimization to be getting the "a" below

        pvalues = _calc_pvalues(junctions, corrected_a, index)
        ecdf = _get_ecdf(pvalues)
        score = _score_ecdf(ecdf)
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

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef float fit_nb(np.ndarray[DTYPE_t, ndim=2] junctionl, str outpath, float nbdisp=0.1, object logger=None) except -1:

    cdef np.ndarray[np.int_t, ndim=1] indices
    cdef np.ndarray[DTYPE_t, ndim=1] mean_junc, std_junc
    cdef np.ndarray[DTYPE_t, ndim=1] ecdf, pvalues

    cdef np.ndarray[DTYPE_t, ndim=1] precision_values = np.array([0.1, 0.01, 0.001], dtype=np.float64)
    cdef float one_over_r0, b, one_over_r, score, precision
    cdef int i

    if logger and outpath:
        logger.debug("NBFit: Plots will be drawn in %s..." % outpath)

    if junctionl.shape[0] < 10:
        logger.warning("Your dataset is not deep enougth to define an apropiate NB factor. The default 0 is given")
        return 0.0
    if junctionl.shape[0] < 5000:
        junctions = junctionl
    else:
        junctions = junctionl[np.random.choice(junctionl.shape[0], 5000)]

    #junctions = masked_less(junctionl, 0.1)
    mean_junc = junctions.mean(axis=1)
    std_junc = junctions.std(axis=1)

    indices = np.zeros(shape=len(junctions), dtype=np.int)
    for i, jj in enumerate(junctions):
        jji = jj.nonzero()
        indices[i] = np.random.choice(jji[0])
    # linear regression, retrieve the a and the b plus
    one_over_r0, b = np.polyfit(mean_junc, std_junc, 1)

    pvalues = _calc_pvalues(junctions, one_over_r0, indices)
    ecdf = _get_ecdf(pvalues)
    score = _score_ecdf(ecdf)
    one_over_r = one_over_r0

    for i, precision in enumerate(precision_values):
        one_over_r, score, ecdf, pvalues = _adjust_fit(one_over_r, junctions, precision, score, outpath,
                                                       index=indices, logger=logger)
        if logger:
            logger.debug("Corrected to %.5f with precision %s. Current score is %.5f" % (one_over_r, precision, score))
        if i + 1 != len(precision_values):
            #go "up" in the scale so we dont miss better solution
            one_over_r += precision - precision_values[i + 1]
            pvalues = _calc_pvalues(junctions, one_over_r, indices)
            ecdf = _get_ecdf(pvalues)
            score = _score_ecdf(ecdf)

    if logger:
        logger.debug("Calculating the nb_r and nb_p with the new fitted function")

    return one_over_r