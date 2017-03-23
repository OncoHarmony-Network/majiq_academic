
import os
import random
import sys
import numpy as np
from scipy.stats import beta
import majiq.src.adjustdelta as majiq_delta
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.io_utils
from majiq.src.plotting import plot_matrix
from majiq.src.constants import *
from majiq.src.sample import sample_from_junctions

from majiq.src.beta_binomial import betabinom
import scipy.misc

"""
Calculate and manipulate PSI and Delta PSI values
"""
BSIZE = 0.025  # TODO To parameters
BINS = np.arange(0, 1, BSIZE)
# The bins for PSI values. With a BSIZE of 0.025, we have 40 BINS
BINS_CENTER = np.arange(0 + BSIZE / 2, 1, BSIZE)
# The center of the previous BINS. This is used to calculate the mean value of each bin.


def divs_from_bootsamples(lsvs_to_work, fitfunc_r, n_replica, pnorm, m_samples, k_positions, discardzeros,
                          trimborder, debug=False, nbins=40, store_bootsamples=True):

    bsize = 1.0 / float(nbins)
    psi_border = np.arange(0, 1.01, bsize)

    div = []#np.zeros(shape=(len(lsvs_to_work), n_replica), dtype=np.float)
    for lsv_idx, quant_lsv in enumerate(lsvs_to_work):
        num_ways = quant_lsv.coverage[0].shape[0]
        post_cdf = np.zeros(shape=(n_replica, num_ways, psi_border.shape[0]), dtype=np.float)
        for rr in xrange(n_replica):

            alpha_prior, beta_prior = get_prior_params(quant_lsv.type, num_ways)
            m_lsv, var_lsv, s_lsv = sample_from_junctions(junction_list=quant_lsv.coverage[rr],
                                                          m=m_samples,
                                                          k=k_positions,
                                                          discardzeros=discardzeros,
                                                          trimborder=trimborder,
                                                          fitted_one_over_r=fitfunc_r[rr],
                                                          debug=debug)

            for jidx in range(num_ways):
                alpha_0 = alpha_prior[jidx]
                beta_0 = beta_prior[jidx]
                for m in xrange(m_samples):

                    sample = s_lsv[jidx][m]
                    notsample = s_lsv[:, m].sum() - sample
                    post_cdf[rr, jidx] += beta.cdf(psi_border, a=sample + alpha_0, b=notsample + beta_0)

                post_cdf[rr, jidx] /= m_samples

        med = np.median(post_cdf, axis=0, keepdims=True)

        psi = np.log(np.diff(post_cdf, axis=2) + MINVAL)
        psi -= scipy.misc.logsumexp(post_cdf, axis=2, keepdims=True)
        med = np.log(np.diff(med, axis=2) + MINVAL)
        med -= scipy.misc.logsumexp(med, axis=2, keepdims=True)
        ldist = np.log(np.abs(np.exp(psi) - np.exp(med)) + MINVAL)
        div.append(np.exp((scipy.misc.logsumexp(ldist * pnorm, axis=2) - np.log(2)) / pnorm).max(axis=1))

    return div


def calc_weights_from_divs(divs, thresh=0.75, alpha=15., logger=None):
    k_T = np.where(divs.max(1) > thresh)[0]
    if k_T.size > 0:
        iter_k_t = (divs[k_T] == divs[k_T].max(1)[:, np.newaxis]).sum(0)
        alpha_param = float(alpha) / iter_k_t.size
        num_N = k_T.size / iter_k_t.size
        beta_param = alpha - alpha_param
        distr = betabinom(k_T.size, alpha_param, beta_param)
        weights = np.clip(distr.sf(iter_k_t) / distr.sf(num_N), 0, 1)
        logger.debug('|K_T|=%d, N=%d, Alpha=%0.03f, Beta=%0.03f' % (k_T.size, num_N, alpha_param, beta_param), logger)
    else:
        weights = np.ones(divs.shape[1], dtype=float)
    return weights


def empirical_delta_psi(list_lsv, lsv_types, lsv_dict1, lsv_summarized1, lsv_dict2, lsv_summarized2, logger=None):
    """Simple PSI calculation without involving a dirichlet prior, coming from reads from junctions
    :type list_lsv: object
    """

    delta_psi = []
    delta_psi_ir = []

    for idx, lsv in enumerate(list_lsv):
        # Assuming that the type is the same in all the replicas and groups
        if lsv_types[lsv].endswith('i'):
            delta_psi_res = delta_psi_ir
        else:
            delta_psi_res = delta_psi

        cov1 = lsv_summarized1[:, lsv_dict1[lsv], 0].mean(axis=0)

        # cov1 = [fp[JUNCTIONS_DATASET_NAME][fp['LSVs/%s' % lsv].attrs['coverage']].sum(axis=1) for fp in group1])
        # cov1 = cov1.mean(axis=0)
        psi1 = np.array([float(cov1[jidx]) / float(np.sum(cov1)) for jidx in range(len(cov1))])
        psi1[np.isnan(psi1)] = 0.5

        cov2 = lsv_summarized2[:, lsv_dict2[lsv], 0].mean(axis=0)
        # cov2 = np.array([fp[JUNCTIONS_DATASET_NAME][fp['LSVs/%s' % lsv].attrs['coverage']].sum(axis=1) for fp in group2])
        # cov2 = cov2.mean(axis=0)
        psi2 = np.array([float(cov2[jidx]) / float(np.sum(cov2)) for jidx in range(len(cov2))])
        psi2[np.isnan(psi2)] = 0.5

        delta_psi_res.append(psi1 - psi2)
        #   if logger: logger.info("Calculating delta PSI for 'best set'...")

    return delta_psi, delta_psi_ir


def __load_default_prior():

    encoding = sys.getfilesystemencoding()
    direc = os.path.dirname(unicode(__file__, encoding))
    def_mat = majiq.src.io_utils.load_bin_file('%s/../data/defaultprior.pickle' % direc)
    return def_mat


def gen_prior_matrix(lsv_dict1, lsv_summarized1, lsv_dict2, lsv_summarized2, lsv_types, output, conf, numbins=20,
                     defaultprior=False, minpercent=-1, logger=None):
    #Start prior matrix
    logger.info("Calculating prior matrix...")
    psi_space = np.linspace(0, 1 - conf.binsize, num=numbins) + conf.binsize / 2
    if defaultprior:
        def_mat = __load_default_prior()
        prior_matrix = [def_mat, def_mat]
        return psi_space, prior_matrix

    logger.debug('Filtering to obtain "best set"...')

    # temp_files = [files[0],files[]]

    filtered_lsv1 = majiq_filter.merge_files_hdf5(lsv_dict1, lsv_summarized1, minnonzero=10, min_reads=20,
                                                  percent=minpercent, logger=logger)
    filtered_lsv2 = majiq_filter.merge_files_hdf5(lsv_dict2, lsv_summarized2, minnonzero=10, min_reads=20,
                                                  percent=minpercent, logger=logger)

    list_of_lsv = list(set(filtered_lsv1).intersection(set(filtered_lsv2)))
    logger.debug("'Best set' is %s events" % len(list_of_lsv))
    best_dpsi, best_dpsi_ir = empirical_delta_psi(list_of_lsv, lsv_types,
                                                  lsv_dict1, lsv_summarized1,
                                                  lsv_dict2, lsv_summarized2, )

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
            mixture_pdf = majiq_delta.adjustdelta_lsv(best_delta_psi, output, plotpath=conf.plotpath,
                                                      title=" ".join(conf.names), numiter=conf.iter,
                                                      breakiter=conf.breakiter, njunc=nj, logger=logger)
            pmat = []
            for i in xrange(numbins):
                pmat.extend(mixture_pdf[numbins - i:(numbins * 2) - i])

            prior_matrix[prior_idx] = np.array(pmat).reshape(numbins, -1)
            if np.isnan(prior_matrix[prior_idx]).any():
                if prior_idx == 1:
                    logger.WARNING("Not enought statistic power to calculate the intron retention specific prior, "
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
                        "prior_matrix_jun_%s" % nj, conf.plotpath)

    return psi_space, prior_matrix


def prob_data_sample_given_psi(sample, all_sample, nbins, alpha_prior, beta_prior):
    bsize = 1.0 / float(nbins)
    psi_border = np.arange(0, 1.01, bsize)
    notsample = all_sample - sample

    bincdf = beta.cdf(psi_border, a=sample + alpha_prior, b=notsample + beta_prior)
    bin_test = bincdf[1:] - bincdf[:-1] + 1e-300

    return bin_test


def combine_for_priormatrix(group1, group2, matched_info, num_exp):
    res_group1 = []
    res_group2 = []

    for lidx, lsv in enumerate(matched_info):
        idx = random.randrange(num_exp[0])
        res_group1.append(group1[lidx][idx])

        idx = random.randrange(num_exp[1])
        res_group2.append(group2[lidx][idx])

    grp1 = [res_group1, matched_info]
    grp2 = [res_group2, matched_info]

    return grp1, grp2


def get_prior_params(lsvinfo, num_ways):
    if 'i' in lsvinfo[2]:
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