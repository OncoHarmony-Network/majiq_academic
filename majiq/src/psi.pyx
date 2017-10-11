import os
import sys
from scipy.stats import beta
import majiq.src.adjustdelta as majiq_delta
import majiq.src.filter as majiq_filter
from majiq.src.plotting import plot_matrix
from majiq.src.constants import *
from majiq.src.beta_binomial import betabinom
import scipy.misc
import numpy as np
cimport numpy as np
import cython
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
@cython.wraparound(True)  # turn off negative index wrapping for entire function
cdef np.ndarray _prob_data_sample_given_psi(float sample, float all_sample, int nbins, float alpha_prior,
                                            float beta_prior):
    cdef float bsize
    cdef np.ndarray psi_border, bincdf, bin_test
    cdef float notsample

    bsize = 1.0 / float(nbins)
    psi_border = np.arange(0, 1.01, bsize)
    notsample = all_sample - sample

    bincdf = beta.cdf(psi_border, a=sample + alpha_prior, b=notsample + beta_prior)
    bin_test = bincdf[1:] - bincdf[:-1] + 1e-300

    return bin_test

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(True)  # turn off negative index wrapping for entire function
cdef np.ndarray _samples_from_psi(np.ndarray post_psi, float mu_psi, float vwindow, int nsamples, int nbins):
    cdef float bsize
    cdef np.ndarray psi_border, outsamps, cdf, pr, mask
    cdef list xs

    bsize = 1.0 / float(nbins)
    psi_border = np.arange(0, 1.01, bsize)
    outsamps = np.random.choice(psi_border[:-1], p=post_psi, size=nsamples)
    outsamps += np.random.rand(nsamples) * bsize
    if vwindow > 0:
        xs = [mu_psi - vwindow, mu_psi + vwindow]
        cdf = np.append([0.], post_psi.cumsum())
        pr = np.interp(xs, psi_border, cdf, left=0., right=1.).dot([-1, 1])
        if pr == 0:
            outsamps = -1.
        elif pr < 1:
            mask = np.random.rand(nsamples) < pr
            outsamps *= mask
            outsamps += mask - 1
    return outsamps

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple _psi_posterior(psi, int p_idx, int m_samples, int num_exp, int num_ways, int nbins,
                    float alpha_0, float beta_0):

    cdef np.ndarray posterior = np.zeros(shape=nbins, dtype=np.float)
    cdef np.ndarray data_given_psi
    cdef list mu_psi_m = []
    cdef float mu_psi
    cdef int m
    cdef np.ndarray alls
    alls = psi.sum(axis=(0,1))
    for m in range(m_samples):
        # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
        junc = psi[:, p_idx, m]
        #all_sample = np.array([psi[xx][yy][m].sum() for xx in range(num_exp) for yy in range(num_ways)])
        mu_psi_m.append(float(junc.sum()) / alls[m].sum())
        data_given_psi = np.log(_prob_data_sample_given_psi(junc.sum(), alls[m].sum(), nbins, alpha_0, beta_0))
        posterior += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))

    mu_psi = np.median(mu_psi_m)
    posterior /=  m_samples

    return mu_psi, posterior

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple _deltapsi_posterior(np.ndarray psi1, np.ndarray psi2, np.ndarray prior_matrix, int p_idx, int m_samples,
                         list num_exp, int num_ways, int nbins, float alpha_0, float beta_0):

    cdef np.ndarray ones_n = np.ones(shape=(1, nbins), dtype=np.float)
    cdef np.ndarray posterior = np.zeros(shape=(nbins, nbins), dtype=np.float)
    cdef np.ndarray post_psi1 = np.zeros(shape=nbins, dtype=np.float)
    cdef np.ndarray post_psi2 = np.zeros(shape=nbins, dtype=np.float)
    cdef np.ndarray A, psi_v1, psi_v2
    cdef np.ndarray data_given_psi1, data_given_psi2
    cdef list mu_psi1_m = []
    cdef list mu_psi2_m = []
    cdef int m

    cdef float junc
    cdef float all_sample
    cdef float mu_psi1, mu_psi2
    cdef np.ndarray alls1, alls2
    print(psi1.shape)
    alls1 = psi1.sum(axis=(0,1))
    alls2 = psi2.sum(axis=(0,1))

    for m in range(m_samples):
        # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))

        junc = psi1[:, p_idx, m].sum()
        #TODO: check this

        # all_sample = np.array([psi1[xx, yy, m].sum() for xx in range(num_exp[0]) for yy in range(num_ways)]).sum()
        # if all_sample - alls1[m] > 0.00001: print (all_sample, psi1[:,:,m].sum(), alls1[m])
        mu_psi1_m.append(float(junc + alpha_0) / (alls1[m] + alpha_0 + beta_0))
        data_given_psi1 = np.log(_prob_data_sample_given_psi(junc, alls1[m], nbins, alpha_0, beta_0))

        psi_v1 = data_given_psi1.reshape(nbins, -1)
        post_psi1 += np.exp(data_given_psi1 - scipy.misc.logsumexp(data_given_psi1))

        junc = psi2[:, p_idx, m].sum()
        # all_sample = np.array([psi2[xx][yy][m] for xx in range(num_exp[1]) for yy in range(num_ways)]).sum()
        # if all_sample - alls2[m] > 0.00001: print (all_sample, psi2[:,:,m].sum(), alls2[m])
        mu_psi2_m.append(float(junc + alpha_0) / (alls2[m] + alpha_0 + beta_0))
        data_given_psi2 = np.log(_prob_data_sample_given_psi(junc, alls2[m], nbins, alpha_0, beta_0))

        post_psi2 += np.exp(data_given_psi2 - scipy.misc.logsumexp(data_given_psi2))
        psi_v2 = data_given_psi2.reshape(-1, nbins)

        A = (psi_v1 * ones_n + psi_v2 * ones_n.T) + np.log(prior_matrix)
        posterior += np.exp(A - scipy.misc.logsumexp(A))


    mu_psi1 = np.median(mu_psi1_m)
    mu_psi2 = np.median(mu_psi2_m)
    posterior /= m_samples
    post_psi1 /= m_samples
    post_psi2 /= m_samples

    return mu_psi1, mu_psi2, posterior, post_psi1, post_psi2

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple _heterogen_posterior(np.ndarray boots, int m_samples, int psi_samples, float vwindow, int num_exp, int nbins,
                          np.ndarray alpha_prior, np.ndarray beta_prior):
    cdef int num_ways =  boots.shape[1]

    mu_psi = np.zeros(shape=(num_exp, num_ways))
    mean_psi = np.zeros(shape=(num_ways, nbins), dtype=np.float)
    samps = np.zeros(shape=(num_exp, num_ways, psi_samples))

    for exp in range(num_exp):
        all_sample = boots[exp].sum(axis=0)
        for p_idx in range(num_ways):
            alpha_0 = alpha_prior[p_idx]
            beta_0 = beta_prior[p_idx]
            post_psi = np.zeros(shape=nbins, dtype=np.float)
            for m in range(m_samples):
                junc = boots[exp, p_idx, m]
                data_given_psi = np.log(_prob_data_sample_given_psi(junc, all_sample[m], nbins, alpha_0, beta_0))
                post_psi += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))
                mu_psi[exp, p_idx] += float(junc + alpha_0) / (all_sample[m] + alpha_0 + beta_0)

            post_psi /= m_samples
            mean_psi[p_idx] += post_psi
            mu_psi[exp, p_idx] /= m_samples

            if psi_samples == 1:
                samps[exp, p_idx, 0] = mu_psi
            else:
                samps[exp, p_idx, :] = _samples_from_psi(post_psi, mu_psi[exp, p_idx], vwindow, psi_samples, nbins)

    mean_psi /= num_exp

    return samps, mean_psi, mu_psi

# External calls
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef tuple psi_posterior(np.ndarray psi, int m, int num_exp, int nbins, str lsv_type):

    cdef list mu_psi = []
    cdef list post_psi = []
    cdef int num_ways = <int> psi.shape[1]
    cdef int p_idx
    cdef np.ndarray alpha_prior, beta_prior
    alpha_prior, beta_prior = _get_prior_params(lsv_type, num_ways)

    for p_idx in range(int(num_ways)):
        mu_psi_m, posterior = _psi_posterior(psi, p_idx, m, num_exp, num_ways, nbins, alpha_prior[p_idx],
                                             beta_prior[p_idx])

        mu_psi.append(mu_psi_m)
        post_psi.append(posterior)
        if num_ways == 2:
            break

    return mu_psi, post_psi

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef tuple deltapsi_posterior(np.ndarray psi1, np.ndarray psi2, np.ndarray prior_matrix, int m, list num_exp, int nbins,
                         str lsv_type):

    cdef list mu_psi1 = []
    cdef list mu_psi2 = []
    cdef list post_matrix = []
    cdef list posterior_psi1 = []
    cdef list posterior_psi2 = []
    cdef int num_ways = <int> psi1.shape[1]
    cdef np.ndarray alpha_prior, beta_prior
    cdef int p_idx, prior_idx
    cdef tuple vals

    alpha_prior, beta_prior = _get_prior_params(lsv_type, num_ways)
    prior_idx = 1 if 'i' in lsv_type else 0
    for p_idx in range(num_ways):
        vals = _deltapsi_posterior(psi1, psi2, prior_matrix[prior_idx], p_idx, m, num_exp, num_ways, nbins,
                                   alpha_prior[p_idx], beta_prior[p_idx])

        mu_psi1.append(vals[0])
        mu_psi2.append(vals[1])
        post_matrix.append(vals[2])
        posterior_psi1.append(vals[3])
        posterior_psi2.append(vals[4])

        if num_ways == 2:
            break

    return post_matrix, posterior_psi1, posterior_psi2, mu_psi1, mu_psi2

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef list heterogen_posterior(list boots, object out_het, int m_samples, int psi_samples, float vwindow, list num_exp,
                          int nbins, str lsv_type):

    cdef np.ndarray alpha_prior, beta_prior
    cdef list samps
    cdef int num_ways =  boots[0].shape[1]

    samps = [None, None]
    alpha_prior, beta_prior = _get_prior_params(lsv_type, num_ways)
    for grp_idx in range(2):
        samps[grp_idx], mean_psi, mu_psi = _heterogen_posterior(boots[grp_idx], m_samples, psi_samples, vwindow,
                                                                num_exp[grp_idx], nbins, alpha_prior, beta_prior)

        out_het.add_group(mu_psi, mean_psi)
        #print(grp_idx, mu_psi)

    return samps

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef list divs_from_bootsamples(list lsvs_to_work, int n_replica, float pnorm, int nbins=40):

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
        alpha_prior, beta_prior = _get_prior_params(quant_lsv.type, num_ways)

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
            p_given_signal = rho.dot(p_given_rep)
            cur_wt[ii, rep] = np.clip(p_given_signal / p_given_rep[rep], 0, 1)
    cur_wt = np.array(cur_wt)
    return cur_wt.squeeze()

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple _empirical_delta_psi(list list_of_lsv, dict lsv_empirical_psi1, dict lsv_empirical_psi2, dict lsv_graphic):
    """
    Simple PSI calculation without involving a dirichlet prior, coming from reads from junctions
    """
    cdef str lsv
    cdef list delta_psi = []
    cdef list delta_psi_ir = []

    for lsv in list_of_lsv:
        # Assuming that the type is the same in all the replicas and groups
        if lsv_graphic[lsv].lsv_type.endswith('i'):
            delta_psi_res = delta_psi_ir
        else:
            delta_psi_res = delta_psi
        delta_psi_res.append(lsv_empirical_psi1[lsv] - lsv_empirical_psi2[lsv])

    return delta_psi, delta_psi_ir


def __load_default_prior():

    encoding = sys.getfilesystemencoding()
    direc = os.path.dirname(__file__)

    fop = open('%s/../data/defaultprior.pickle' % direc, 'rb')
    fast_pickler = pickle.Unpickler(fop)
    data = fast_pickler.load()
    fop.close()

    return data


def gen_prior_matrix(dict lsv_dict_graph, dict lsv_empirical_psi1, dict lsv_empirical_psi2, str output, list names,
                     float breakiter, str plotpath, int iter, float binsize, int numbins=20, bint defaultprior=False,
                     int minpercent=-1, object logger=None):

    cdef np.ndarray psi_space, def_mat, lsv, mixture_pdf
    cdef list prior_matrix, list_of_lsv, best_dpsi, best_dpsi_ir, njun_prior, pmat
    cdef int prior_idx, nj


    #Start prior matrix
    logger.info("Calculating prior matrix...")
    psi_space = np.linspace(0, 1 - binsize, num=numbins) + binsize / 2
    if defaultprior:
        def_mat = __load_default_prior()
        prior_matrix = [def_mat, def_mat]
        return psi_space, prior_matrix

    logger.debug('Filtering to obtain "best set"...')

    # temp_files = [files[0],files[]]

    list_of_lsv = list(set(lsv_empirical_psi1.keys()).intersection(set(lsv_empirical_psi2.keys())))
    logger.debug("'Best set' is %s events" % len(list_of_lsv))
    best_dpsi, best_dpsi_ir = _empirical_delta_psi(list_of_lsv, lsv_empirical_psi1, lsv_empirical_psi2, lsv_dict_graph)

    prior_matrix = [[], []]
    for prior_idx, best_delta_psi in enumerate((best_dpsi, best_dpsi_ir)):
        njun_prior = [[]]

        for lsv in best_delta_psi:
            if lsv.shape[0] != 2:
                continue
            njun_prior[0].append(lsv[0])

        for nj in range(len(njun_prior)):

            best_delta_psi = np.array(njun_prior[nj])
            if len(best_delta_psi) == 0:
                if prior_idx == 0:
                    prior_matrix[prior_idx] = __load_default_prior()
                else:
                    prior_matrix[prior_idx] = prior_matrix[0]
                continue

            logger.debug("Parametrizing 'best set'...%s", prior_idx)
            mixture_pdf = majiq_delta.adjustdelta_lsv(best_delta_psi, output, plotpath=plotpath,
                                                      title=" ".join(names), numiter=iter,
                                                      breakiter=breakiter, njunc=nj, logger=logger)
            pmat = []
            for i in range(numbins):
                pmat.extend(mixture_pdf[numbins - i:(numbins * 2) - i])

            prior_matrix[prior_idx] = np.array(pmat).reshape(numbins, -1)
            if np.isnan(prior_matrix[prior_idx]).any():
                if prior_idx == 1:
                    logger.warning("Not enought statistic power to calculate the intron retention specific prior, "
                                   "in that case we will use the global prior")
                    prior_matrix[prior_idx] = prior_matrix[0]
                else:
                    raise ValueError(" The input data does not have enought statistic power in order to calculate "
                                     "the prior. Check if the input is correct or use the --default_prior option in "
                                     " order to use a precomputed prior")
            else:
                prior_matrix[prior_idx] /= sum(prior_matrix[prior_idx])
                # renormalize so it sums 1

            plot_matrix(prior_matrix[prior_idx], "Prior Matrix , version %s" % prior_idx,
                        "prior_matrix_jun_%s" % nj, plotpath)

    return psi_space, prior_matrix


