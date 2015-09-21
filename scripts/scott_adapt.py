import sys
import scipy.misc
import numpy as np
import pickle
from scipy.stats import beta


def __prob_data_sample_given_psi(sample, all_sample, nbins, alpha_prior, beta_prior):
    bsize = 1.0 / float(nbins)
    psi_border = np.arange(0, 1.01, bsize)
    notsample = all_sample - sample
    bincdf = [beta.cdf(xx, a=sample + alpha_prior, b=notsample + beta_prior) for xx in psi_border]

    bin_test = []
    for x in xrange(nbins):
        val = bincdf[x + 1] - bincdf[x]
        bin_test.append(val)

    bin_test = np.array(bin_test) + 1e-300

    return bin_test


def __collapse_matrix(matrix):
    """Collapse the diagonals probabilities in 1-D and return them"""
    collapse = []
    matrix_corner = matrix.shape[0]
    for i in xrange(-matrix_corner+1, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)


def __get_prior_params(lsvinfo, num_ways):
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


def calc_psi(info, lsv_samples, num_exp, weights,
             nbins=40, msamples=100):

    """
    :param info: list of metainfo for lsv, the order is the same than lsv_samples
    :param lsv_samples: list of lsv samples. each lsv junction has msamples samples
    :param num_exp: number of replicas to quantify
    :param nbins: number of bins for the posterior marginalization
    :param msamples: number of samples that we get from the bootstraping
    :return: tuple with two lsv list in the same order (posterior distribution, metainformation)
             This posterior distribution is one matrix per junction in the lsv that specifies the distribution. In
             order to collapse the matrix to a nbins histogram, you should execute the function collapse_matrix
    """

    post_psi = []
    new_info = []
    for lidx, lsv_info in enumerate(info):
        num_ways = len(lsv_samples[lidx, 0])
        if lidx % 50 == 0:
            print "Event %d ..." % lidx
            sys.stdout.flush()

        alpha_prior, beta_prior = __get_prior_params(lsv_info, num_ways)

        psi = lsv_samples[lidx, :]
        post_psi.append([])
        new_info.append(lsv_info)
        for p_idx in xrange(int(num_ways)):
            posterior = np.zeros(shape=nbins, dtype=np.float)

            for m in xrange(msamples):
                # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                junc = [psi[xx][p_idx][m] for xx in xrange(num_exp)]
                junc = np.array(junc)
                all_sample = [psi[xx][yy][m].sum() for xx in xrange(num_exp) for yy in xrange(num_ways)]
                all_sample = np.array(all_sample)
                data_given_psi = np.log(__prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                     alpha_prior[p_idx], beta_prior[p_idx]))

                # normalizing
                posterior += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))

            post_psi[-1].append(posterior / msamples)
            if num_ways == 2:
                break

    return post_psi, new_info


def deltapsi(info, lsv_samples1, lsv_samples2, num_exp, weights,
             prior_path='./data/default_prior.pickle', nbins=40, msamples=100):
    """

    :param info: list of metainfo for lsv, the order is the same than lsv_samples
    :param lsv_samples1: list of lsv samples for condition 1. each lsv junction has msamples sample.
    :param lsv_samples2: list of lsv samples for condition 2. each lsv junction has msamples sample.
    :param num_exp: number of replicas to quantify
    :param prior_path: path to the pickle file with the prior distribution
    :param nbins: number of bins for the posterior marginalization
    :param msamples: number of samples that we get from the bootstraping
    :return: tuple with 4 lsv list in the same order (posterior distribution, metainformation, psi1  and psi2)
             +posterior distribution: one matrix per junction in the lsv that specifies the distribution. In
             order to collapse the matrix to a nbins histogram, you should execute the function collapse_matrix
             +metainformation: the information of the lsv, similar to info
             +psi1: the psi posterior distribution for condition 1, the same as result of calcpsi
             +psi1: the psi posterior distribution for condition 2, the same as result of calcpsi

    """

    post_matrix = []
    new_info = []
    ones_n = np.ones(shape=(1, nbins), dtype=np.float)

    prior_matrix = pickle.load(open(prior_path))

    posterior_psi1 = []
    posterior_psi2 = []

    for lidx, lsv_info in enumerate(info):
        num_ways = len(lsv_samples1[lidx][0])
        if lidx % 50 == 0:
            print "Event %d ..." % lidx,
            sys.stdout.flush()

        alpha_prior, beta_prior = __get_prior_params(lsv_info, num_ways)
        if 'i' in lsv_info[2]:
            prior_idx = 1
        else:
            prior_idx = 0
        post_matrix.append([])
        new_info.append(lsv_info)

        posterior_psi1.append([])
        posterior_psi2.append([])
        psi1 = lsv_samples1[lidx, :]
        psi2 = lsv_samples2[lidx, :]

        for p_idx in xrange(num_ways):

            posterior = np.zeros(shape=(nbins, nbins), dtype=np.float)
            post_psi1 = np.zeros(shape=nbins, dtype=np.float)
            post_psi2 = np.zeros(shape=nbins, dtype=np.float)
            for m in xrange(msamples):
                # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                junc = [psi1[xx][p_idx][m] for xx in xrange(num_exp[0])]
                junc = np.array(junc)
                all_sample = np.array([psi1[xx][yy][m].sum() for xx in xrange(num_exp[0]) for yy in xrange(num_ways)])
                data_given_psi1 = np.log(__prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                    alpha_prior[p_idx], beta_prior[p_idx]))

                psi_v1 = data_given_psi1.reshape(nbins, -1)
                post_psi1 += np.exp(data_given_psi1 - scipy.misc.logsumexp(data_given_psi1))

                junc = [psi2[xx][p_idx][m] for xx in xrange(num_exp[1])]
                junc = np.array(junc)
                all_sample = np.array([psi2[xx][yy][m].sum() for xx in xrange(num_exp[1]) for yy in xrange(num_ways)])
                data_given_psi2 = np.log(__prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                    alpha_prior[p_idx], beta_prior[p_idx]))
                post_psi2 += np.exp(data_given_psi2 - scipy.misc.logsumexp(data_given_psi2))
                psi_v2 = data_given_psi2.reshape(-1, nbins)

                A = (psi_v1 * ones_n + psi_v2 * ones_n.T) + np.log(prior_matrix[prior_idx])
                posterior += np.exp(A - scipy.misc.logsumexp(A))

            post_matrix[-1].append(posterior / msamples)
            posterior_psi1[-1].append(post_psi1 / msamples)
            posterior_psi2[-1].append(post_psi2 / msamples)
            if num_ways == 2:
                break

    return post_matrix, new_info, posterior_psi1, posterior_psi2