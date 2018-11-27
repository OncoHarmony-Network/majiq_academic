import os
import sys
from scipy.stats import beta
import majiq.src.adjustdelta as majiq_delta
from majiq.src.internals.psi cimport psi_distr_t
from majiq.src.constants import *
from majiq.src.beta_binomial import betabinom
from libcpp.vector cimport vector
from libcpp.pair cimport pair


import scipy.misc
import numpy as np
cimport numpy as np
import cython
from libcpp.string cimport string
import pickle

"""
Calculate and manipulate PSI and Delta PSI values
"""

# Internal functions
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(True)  # turn off negative index wrapping for entire function
cdef tuple _get_prior_params(str lsv_type, int num_ways):
    cdef float alpha, bta
    cdef np.ndarray alpha_prior, beta_prior

    if 'i' in lsv_type:
        alpha = 1.0 / (num_ways - 1)
        alpha *= float(num_ways) / (num_ways + 1)
        alpha_prior = np.array([alpha] * num_ways)
        alpha_prior[-1] = 1.0 / (num_ways + 1)
        beta_prior = 1 - alpha_prior

    else:
        alpha = 1.0 / num_ways
        bta = float(num_ways - 1.0) / num_ways
        alpha_prior = np.array([alpha] * num_ways)
        beta_prior = np.array([bta] * num_ways)

    return alpha_prior, beta_prior


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef list divs_from_bootsamples(list lsvs_to_work, int n_replica, object lsv_types, float pnorm, int nbins=40):

    cdef float bsize, alpha_0, beta_0
    cdef np.ndarray psi_border, psi, med, post_cdf, ldist
    cdef list div = []
    cdef int lsv_idx, rr, num_ways, jidx
    cdef np.ndarray alpha_prior, beta_prior, s_lsv
    cdef float sample, notsample
    cdef int m_samples = 100


    bsize = 1.0 / float(nbins)
    psi_border = np.arange(0, 1.01, bsize)

    for lsv_idx, quant_lsv in enumerate(lsvs_to_work):
        num_ways = quant_lsv.coverage[0].shape[0]
        post_cdf = np.zeros(shape=(n_replica, num_ways, psi_border.shape[0]), dtype=np.float)
        #print(quant_lsv.id, lsv_types[quant_lsv.id])
        alpha_prior, beta_prior = _get_prior_params(lsv_types[quant_lsv.id], num_ways)

        for rr, s_lsv in enumerate(quant_lsv.coverage):
            m_samples = s_lsv.shape[1]
            for jidx in range(num_ways):
                alpha_0 = alpha_prior[jidx]
                beta_0 = beta_prior[jidx]
                for m in range(m_samples):
                    sample = s_lsv[jidx][m]
                    notsample = s_lsv[:, m].sum() - sample
                    post_cdf[rr, jidx] += beta.cdf(psi_border, a=sample + alpha_0, b=notsample + beta_0)

                post_cdf[rr, jidx] /= m_samples

        med = np.median(post_cdf, axis=0, keepdims=True)

        psi = np.log(np.diff(post_cdf, axis=2) + MINVAL)
        psi -= scipy.misc.logsumexp(psi, axis=2, keepdims=True)
        med = np.log(np.diff(med, axis=2) + MINVAL)
        med -= scipy.misc.logsumexp(med, axis=2, keepdims=True)
        ldist = np.log(np.abs(np.exp(psi) - np.exp(med)) + MINVAL)
        div.append(np.exp((scipy.misc.logsumexp(ldist * pnorm, axis=2) - np.log(2)) / pnorm).max(axis=1))

    return div

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def calc_rho_from_divs(np.ndarray divs, float thresh=0.75, float alpha=15., int nreps=1, object logger=None):

    K_T, bad_reps = np.where(divs > thresh)
    K_T = np.unique(K_T)
    if logger:
        logger.debug('%d' % K_T.size)

    if K_T.size > 0:
        K_t = np.bincount(bad_reps, minlength=nreps)
        if logger:
            logger.debug('%s' % (', '.join(['%d' % xx for xx in K_t])))
        Alpha = float(alpha) / K_t.size
        N = K_T.size / K_t.size
        Beta = alpha - Alpha
        distr = betabinom(K_T.size, Alpha, Beta)
        rho = np.clip(distr.sf(K_t) / distr.sf(N), 0, 1)
        if logger:
            logger.debug('|K_T|=%d, N=%d, Alpha=%0.03f, Beta=%0.03f' % (K_T.size, N,
                                                                 Alpha, Beta))
    else:
        rho = np.ones(nreps)
    return rho

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def calc_local_weights(divs, rho, local):
    dshp = divs.shape
    cur_wt = np.zeros(shape=dshp)
    for ii in range(dshp[0]):
        lsv = divs[ii, :]
        for rep, dist in enumerate(lsv):
            p_given_rep = ((divs >= (dist - local)) & (divs <= (dist + local))).mean(axis=0)
            p_given_signal = rho.dot(p_given_rep) / rho.sum()
            cur_wt[ii, rep] = np.clip(p_given_signal / p_given_rep[rep], 0, 1)
    cur_wt = np.array(cur_wt)
    return cur_wt.squeeze()

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple _empirical_delta_psi(list list_of_lsv, dict lsv_empirical_psi1, dict lsv_empirical_psi2, object lsv_type):
    """
    Simple PSI calculation without involving a dirichlet prior, coming from reads from junctions
    """
    cdef string lsv
    cdef list delta_psi = []
    cdef list delta_psi_ir = []

    for lsv in list_of_lsv:
        if lsv_type[lsv][1] > 2 : continue
        # Assuming that the type is the same in all the replicas and groups
        if lsv_type[lsv][0].endswith(b'i'):
            delta_psi_res = delta_psi_ir
        else:
            delta_psi_res = delta_psi
        delta_psi_res.append(lsv_empirical_psi1[lsv][0] - lsv_empirical_psi2[lsv][0])
    return np.array(delta_psi, dtype=np.float32), np.array(delta_psi_ir, dtype=np.float32)
    # return delta_psi, delta_psi_ir


cpdef void __load_default_prior(vector[vector[psi_distr_t]]& prior_matrix):

    cdef int numbins = prior_matrix[0].size()
    cdef int xx, yy

    encoding = sys.getfilesystemencoding()
    direc = os.path.dirname(__file__)

    fop = open('%s/../data/defaultprior.pickle' % direc, 'rb')
    fast_pickler = pickle.Unpickler(fop)
    data = fast_pickler.load().astype(np.float32)
    fop.close()

    for xx in range(numbins):
        for yy in range(numbins):
            prior_matrix[0][xx][yy] = np.log(data[xx][yy])
            prior_matrix[1][xx][yy] = np.log(data[xx][yy])

    return

#
# cdef print_prior(vector[vector[psi_distr_t]] matrix, int nbins):
#
#     sys.stderr('##MATRIX')
#     for xx in range(nbins):
#         for yy in range(nbins):
#            sys.stderr.write("%s, " % matrix[0][xx][yy])
#         sys.stderr.write("\n")
#
# import sys
# def eprint(*args, **kwargs):
#     print(*args, file=sys.stderr, **kwargs)

cpdef vector[vector[psi_distr_t]] gen_prior_matrix(object lsv_type, dict lsv_empirical_psi1,
                                                                      dict lsv_empirical_psi2, str output, list names,
                                                                      str plotpath, int iter, float binsize,
                                                                      int numbins=20, bint defaultprior=False,
                                                                      int minpercent=-1, object logger=None):

    cdef np.ndarray[np.float32_t] mixture_pdf
    cdef list_of_lsv, njun_prior, pmat
    cdef int prior_idx
    cdef np.ndarray[np.float32_t, ndim=1] best_deltap, best_dpsi, best_dpsi_ir
    cdef np.ndarray[np.float32_t, ndim=1] best_delta_psi
    cdef np.ndarray[np.float32_t, ndim=3] np_pmatrix = np.zeros(shape=(3, numbins, numbins), dtype=np.float32)
    cdef vector[vector[psi_distr_t]] prior_matrix = vector[vector[psi_distr_t]](2, vector[psi_distr_t](numbins, psi_distr_t(numbins, 0)))

    #Start prior matrix
    logger.info("Calculating prior matrix...")
    if defaultprior:
        __load_default_prior(prior_matrix)
        return prior_matrix

    logger.debug('Filtering to obtain "best set"...')


    list_of_lsv = list(set(lsv_empirical_psi1.keys()).intersection(set(lsv_empirical_psi2.keys())))
    logger.debug("'Best set' is %s events" % len(list_of_lsv))
    best_dpsi, best_dpsi_ir = _empirical_delta_psi(list_of_lsv, lsv_empirical_psi1, lsv_empirical_psi2, lsv_type)

    for prior_idx, best_delta_psi in enumerate((best_dpsi, best_dpsi_ir)):
        if len(best_delta_psi) == 0:
            if prior_idx == 0:
                __load_default_prior(prior_matrix)
            else:
                prior_matrix[1] = prior_matrix[0]
            return prior_matrix

        logger.debug("Parametrizing 'best set'...%s", prior_idx)
        mixture_pdf = majiq_delta.adjustdelta(best_delta_psi, iter, output, plotpath=plotpath,
                                                  title=" ".join(names), njunc=2, logger=logger)
        pmat = []
        for i in range(numbins):
            pmat.extend(mixture_pdf[numbins - i:(numbins * 2) - i])

        np_pmatrix[prior_idx] = np.array(pmat, dtype=np.float32).reshape(numbins, -1)
        if np.isnan(np_pmatrix[prior_idx]).any():
            if prior_idx == 1:
                logger.warning("Not enought statistic power to calculate the intron retention specific prior, "
                               "in that case we will use the global prior")
                np_pmatrix[prior_idx] = np_pmatrix[0]
            else:
                raise ValueError(" The input data does not have enought statistic power in order to calculate "
                                 "the prior. Check if the input is correct or use the --default-prior option in "
                                 " order to use a precomputed prior")
        else:
            np_pmatrix[prior_idx] /= sum(np_pmatrix[prior_idx])

            # renormalize so it sums 1

        # plot_matrix(prior_matrix[prior_idx], "Prior Matrix , version %s" % prior_idx,
        #             "prior_matrix_jun_%s" % nj, plotpath)

    for xx in range(numbins):
        for yy in range(numbins):
            prior_matrix[0][xx][yy] = np.log(np_pmatrix[0, xx, yy])
            prior_matrix[1][xx][yy] = np.log(np_pmatrix[1, xx, yy])
    return prior_matrix


