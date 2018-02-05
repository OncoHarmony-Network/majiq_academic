import numpy as np
cimport numpy as np
from scipy.stats import beta
from scipy.misc import logsumexp

cdef float PSEUDO = 1e-300
ctypedef np.float64_t DTYPE_t

cdef np.ndarray[DTYPE_t, ndim=1] _calc_mixture_pdf(np.ndarray[DTYPE_t, ndim=2] beta_param,
                                                   np.ndarray[DTYPE_t, ndim=2] pmix,
                                                   np.ndarray[DTYPE_t, ndim=1] x_pos):

    cdef np.ndarray[DTYPE_t, ndim=1] mixture_pdf = np.zeros(shape=x_pos.shape[0])
    cdef np.ndarray[DTYPE_t, ndim=1] bincdf
    cdef tuple bt
    cdef int ii

    for ii, bt in enumerate(beta_param):

        bincdf = beta.cdf(x_pos, a=bt[0], b=bt[1])
        mixture_pdf = mixture_pdf + (bincdf[1:] - bincdf[:-1] + 1e-300) * pmix[ii]

    return mixture_pdf


def calculate_beta_params(mean, vari, sample_size):
    b = (1 - mean) * sample_size
    a = mean * sample_size

    return [a, b]


cpdef np.ndarray[DTYPE_t, ndim=1] adjustdelta_lsv(np.ndarray[DTYPE_t, ndim=2] deltapsi, int num_iter, str output,
                                                  str plotpath=None, str title=None,int  njunc=1, object logger=False):

    cdef np.ndarray[DTYPE_t, ndim=2] D = np.zeros(shape=(79, 2), dtype=np.float)
    cdef np.ndarray[DTYPE_t, ndim=1] xpos = np.arange(-1 + (0.025 / 2), 1, 0.025)
    cdef np.ndarray[DTYPE_t, ndim=1] x_pos0 = np.arange(0, 1, 0.025)
    cdef np.ndarray[DTYPE_t, ndim=1] p_mixture
    cdef np.ndarray[DTYPE_t, ndim=2] beta_params
    cdef np.ndarray[DTYPE_t, ndim=1] z_mixture_pdf

    cdef int idx, zero_idx
    cdef float ii, total, num_spike, num_lowcenter, spike
    cdef list uniform = [1, 1]
    cdef list labels = ['Uniform', 'center', 'spike']


    for idx, ii in enumerate(xpos[:-1]):
        D[idx, 0] = round((ii + xpos[idx + 1]) / 2, 5)
        if D[idx, 0] == 0:
            zero_idx = idx

    for ppv in deltapsi:
        for idx, ii in enumerate(xpos[:-1]):
            if ii <= ppv < xpos[idx + 1]:
                D[idx, 1] += 1
                break

    total = D[:, 1].sum()
    num_spike = D[zero_idx, 1]
    num_lowcenter = D[zero_idx - 3:zero_idx + 3, 1].sum() - num_spike

    p_mixture = np.zeros(shape=3, dtype=np.float)

    spike = calculate_beta_params(0.5, 0, num_spike)
    p_mixture[2] = num_spike / total

    center = calculate_beta_params(0.5, 0, num_lowcenter)
    p_mixture[1] = num_lowcenter / total
    p_mixture[0] = 1 - ((num_lowcenter + num_spike) / total)
    beta_params = np.array([uniform, center, spike])

    _em_beta_mix(D, p_mixture, beta_params, num_iter, min_ratio=1e-5, logger=logger, plotpath=plotpath, nj=njunc,
                 labels=labels)

    z_mixture_pdf = _calc_mixture_pdf(beta_params, p_mixture, x_pos0)
    return z_mixture_pdf


cdef tuple _loglikelihood(np.ndarray[DTYPE_t, ndim=2] D, np.ndarray[DTYPE_t, ndim=2] beta_mix,
                    float logp_mix, object logger=False):

    ''' logp_DgK = log P (D | model K ) for each data point without the weight '''
    cdef int N = D.shape[0]
    cdef int K = beta_mix.shape[0]
    cdef int k
    cdef np.ndarray[DTYPE_t, ndim=2] logp_DgK = np.zeros(shape=(N, K), dtype=np.float)
    cdef np.ndarray[DTYPE_t, ndim=2] logp_D
    cdef np.ndarray[DTYPE_t, ndim=1] logp_Dsum = np.zeros(shape=N, dtype=np.float)
    cdef np.ndarray[DTYPE_t, ndim=1] dm
    cdef np.ndarray zrow, no_zrow
    cdef float ll_sum

    for k in range(K):
        logp_DgK[:, k] = np.log(beta.pdf(D[:, 0], beta_mix[k, 0], beta_mix[k, 1]) + PSEUDO)

    logp_D = logp_DgK + logp_mix * np.ones(shape=(N, 1), dtype=np.float)

    dm = np.sum(logp_DgK, axis=1)
    zrow = dm.astype(np.bool)
    no_zrow = np.logical_not(zrow)

    logp_Dsum[no_zrow] = logsumexp(logp_D[no_zrow, :], axis=1)
    logp_Dsum[zrow] = logsumexp(logp_D[zrow, :-1], axis=1)

    ll_sum = np.sum(logp_Dsum * D[:, 1], axis=0)

    return logp_D, logp_Dsum, ll_sum, zrow


cdef _em_beta_mix(np.ndarray[DTYPE_t, ndim=2] D, np.ndarray[DTYPE_t, ndim=1] pmix,
                  np.ndarray[DTYPE_t, ndim=2] beta_mix, num_iter, min_ratio=1e-5,
                  plotpath=None, nj=0, labels=None, logger=False):

    cdef np.ndarray[DTYPE_t, ndim=2] D0 = D.copy()
    cdef np.ndarray[DTYPE_t, ndim=2] logp_D
    cdef np.ndarray[DTYPE_t, ndim=1] logp_Dsum
    cdef np.ndarray zrow
    cdef np.ndarray[DTYPE_t, ndim=1] logp_mix = np.log(pmix)

    cdef np.ndarray[DTYPE_t, ndim=1] new_pmix
    cdef np.ndarray[DTYPE_t, ndim=2] new_beta_mix
    cdef np.ndarray[DTYPE_t, ndim=1] avgxPerK, avgx2PerK, varxPerK
    cdef int N = D.shape[0]
    cdef int K = beta_mix.shape[0]
    cdef float ll_sum, ll_sum_old

    cdef int c = 1
    cdef int a = -1

    if min(D[:, 0]) < 0.0:
        D[:, 0] = (D[:, 0] + 1) / (c - a)

    logp_D, logp_Dsum, ll_sum, zrow = _loglikelihood(D, beta_mix, logp_mix)


    if logger:
        logger.debug("[NJ:%s] Initial Log_Likelihood %.3f \n" % (nj, ll_sum))

    ones_1k = np.ones(shape=(1, K), dtype=np.float)
    for mm in range(num_iter):
        new_beta_mix = np.zeros(shape=(beta_mix.shape[0], beta_mix.shape[1]), dtype=np.float)

        ''' E STEP: '''
        p_KgD = np.exp(logp_D - (logp_Dsum * ones_1k.T).T)
        p_KgD[zrow, K - 1] = 0
        #        pdb.set_trace()

        avgxPerK = np.sum(p_KgD * (D[:, 0] * ones_1k.T).T * (D[:, 1] * ones_1k.T).T, axis=0) / np.sum(
            p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        avgx2PerK = np.sum(p_KgD * np.square((D[:, 0] * ones_1k.T).T) * (D[:, 1] * ones_1k.T).T, axis=0) / np.sum(
            p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        varxPerK = avgx2PerK - (np.square(avgxPerK))

        new_beta_mix[:, 0] = avgxPerK * (((avgxPerK * (1 - avgxPerK)) / varxPerK) - 1)
        new_beta_mix[:, 1] = (1 - avgxPerK) * (((avgxPerK * (1 - avgxPerK)) / varxPerK) - 1)

        new_pmix = np.sum(p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        new_pmix = new_pmix / np.sum(new_pmix, axis=0)

        ll_sum_old = ll_sum
        logp_D, logp_Dsum, ll_sum, zrow = _loglikelihood(D, new_beta_mix, np.log(new_pmix))

        if logger: logger.debug("[NJ:%s] EM Iteration %d:\t ll_sum: %.3f\n" % (nj, mm, ll_sum))

        if ll_sum < ll_sum_old:
            if logger:
                logger.debug("Log_Likelihood DECREASE new %d old %d - Aborting ....\n" % (ll_sum, ll_sum_old))
            break

        pmix = new_pmix.copy()
        beta_mix = new_beta_mix.copy()
        logp_mix = np.log(pmix)

        if np.exp(ll_sum - ll_sum_old) < (1.0 + min_ratio):
            if logger:
                logger.debug("Ratio = %.3f < 1+R(%.3f) - Aborting ... \n" % (ll_sum - ll_sum_old, min_ratio))
            break

