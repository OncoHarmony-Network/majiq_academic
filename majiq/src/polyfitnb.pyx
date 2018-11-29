import numpy as np
cimport numpy as np
from scipy.stats import nbinom, poisson
import cython
from libcpp.vector cimport vector

ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray[np.float64_t, ndim=1] _get_ecdf(np.ndarray[np.float64_t, ndim=1] pvalues):
    cdef int nbins
    cdef np.ndarray hist, bin_edges

    # print sorted(pvalues)
    nbins = max(min(10, len(pvalues)), int(len(pvalues) / 10))
    hist, bin_edges = np.histogram(pvalues, range=[0, 1], bins=nbins, density=True)
    return (np.cumsum(hist) / len(bin_edges))


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef float _score_ecdf(np.ndarray[DTYPE_t, ndim=1] ecdf):
    """
    Give a score to a ecdf calculating the deviation from the 45 degree line
    """
    return sum(abs(np.linspace(0, 1, num=len(ecdf)) - ecdf))


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray _calc_pvalues(np.ndarray[np.float64_t, ndim=2] junctions, float one_over_r,
                              np.ndarray[np.int_t, ndim=1] indices_list):
    cdef np.ndarray[np.float64_t, ndim=1] pvalues,
    cdef int njuncs, idx
    cdef np.ndarray[np.float64_t, ndim=1] junc_idxs, mu
    cdef np.ndarray[np.float64_t, ndim=2] junc_fltr
    cdef np.ndarray[np.int_t, ndim=1] vals
    # cdef np.float32_t r, p

    vals = np.count_nonzero(junctions, axis=1)
    junc_fltr = junctions[vals > 1]
    junc_idxs = np.array([xx[indices_list[idx]] for idx, xx in enumerate(junc_fltr)])
    vals = vals[vals > 1] - 1
    njuncs = vals.shape[0]
    mu = ((junc_fltr.sum(axis=1, dtype=np.float) - junc_idxs) / vals)
    if one_over_r > 0:
        r = 1 / one_over_r
        p = r/ (mu +r)
        pvalues = nbinom.cdf(junc_idxs, r, p)
    else:
        pvalues = poisson.cdf(junc_idxs, mu)
    return 1 - pvalues



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple _adjust_fit(float starting_a, np.ndarray[DTYPE_t, ndim=2] junctions, float precision,
                       float previous_score, np.ndarray[np.int_t, ndim=1] index, object logger=None):

    cdef int previous_a = -1
    cdef int idx = 0
    cdef np.ndarray steps = np.arange(starting_a, 0, - precision)
    cdef np.ndarray[DTYPE_t, ndim=1] pvalues, previous_pvalues, ecdf, previous_ecdf
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


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
# cpdef float fit_nb(vector[np.float32_t *] junctionl, int total_juncs, int eff_length, float nbdisp=0.1, object logger=None) except -1:
cdef float fit_nb(vector[np.float32_t *] junctionl, int total_juncs, int eff_length, float nbdisp, object logger):
    cdef np.ndarray[np.int64_t, ndim=1] indices
    cdef np.ndarray[DTYPE_t, ndim=1] mean_junc, std_junc
    cdef np.ndarray[DTYPE_t, ndim=1] ecdf, pvalues
    cdef np.ndarray[DTYPE_t, ndim=2] junctions


    cdef np.ndarray[DTYPE_t, ndim=1] precision_values = np.array([0.1, 0.01, 0.001], dtype=np.float)
    cdef float one_over_r0, b, one_over_r, score, precision
    cdef int i, j, ix
    cdef DTYPE_t c
    cdef int array_size
    cdef vector[DTYPE_t *] jlist


    if total_juncs < 10:
        logger.warning("Your dataset is not deep enougth to define an apropiate NB factor. The default 0 is given")
        return 0.0

    junctions = np.zeros(shape=(total_juncs, eff_length), dtype=np.float)

    for i in range(total_juncs):
        for j in range(eff_length):
            junctions[i,j] = junctionl[i][j]

    junctions = junctions[junctions.sum(axis=1)>=10]
    if junctions.shape[0] < 10:
        logger.warning("Your dataset is not deep enougth to define an apropiate NB factor. The default 0 is given")
        return 0.0
    array_size = junctions.shape[0] if junctions.shape[0] < 5000 else 5000
    junctions = junctions[np.random.choice(junctions.shape[0], array_size, replace=False)]

    # while ix < total_juncs and i < array_size:
    #     ix +=1
    #     c = 0
    #     for j in range(eff_length):
    #         junctions[i, j] = junctionl[i][j]
    #         c += junctionl[i][j]
    #     if c == 0: continue
    #     i += 1
    #
    # if junctions.sum() == 0.0 :
    #     return 0.0

    junctions[junctions == 0] = np.nan
    mean_junc = np.nanmean(junctions, axis=1)
    std_junc = np.nanstd(junctions, axis=1)
    one_over_r0, b = np.polyfit(mean_junc, std_junc, 1)
    junctions = np.nan_to_num(junctions)
    #[junctions == np.nan] = 0
    indices = np.zeros(shape=len(junctions), dtype=np.int)
    for i, jj in enumerate(junctions):
        jji =jj.nonzero()
        indices[i] = np.random.choice(jji[0])
    # linear regression, retrieve the a and the b plus
    pvalues = _calc_pvalues(junctions, one_over_r0, indices)
    ecdf = _get_ecdf(pvalues)
    score = _score_ecdf(ecdf)
    one_over_r = one_over_r0
    for i, precision in enumerate(precision_values):
        logger.debug(' [Step %s] 1/r: %s' %(i, one_over_r))
        one_over_r, score, ecdf, pvalues = _adjust_fit(one_over_r, junctions, precision, score, index=indices,
                                                       logger=logger)
        if logger:
            logger.debug(" [Step %s] Corrected to %.5f with precision %s. Current score is %.5f" % (i, one_over_r,
                                                                                                    precision, score))
        if i + 1 != len(precision_values):
            #go "up" in the scale so we dont miss better solution
            one_over_r += precision - precision_values[i + 1]
            pvalues = _calc_pvalues(junctions, one_over_r, indices)
            ecdf = _get_ecdf(pvalues)
            score = _score_ecdf(ecdf)

        if logger:
            logger.debug(" [Step %s]  Calculating the nb_r and nb_p with the new fitted function" % i)

    logger.debug('Chosen 1/r NB factor %s' % one_over_r)
    return one_over_r
