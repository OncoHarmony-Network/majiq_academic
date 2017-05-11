import argparse
import pickle
import sys

import numpy as np
import scipy.misc

import majiq.src.filter as majiq_filter
import majiq.src.psi as majiq_psi
from voila.io_voila import VoilaInput
from voila.vlsv import VoilaLsv

BINSIZE = 0.025
NUMBINS = 20

def _dump_lsvs_voila(pickle_path, posterior_matrix, lsvs_info, meta_info, psi_list1=None, psi_list2=None):
    """Create VoilaLSVs objects readable by voila."""
    vlsvs = []
    psi1, psi2 = None, None
    for ii, bins in enumerate(posterior_matrix):
        lsv_graphic = lsvs_info[ii][-1]
        if psi_list1:
            psi1, psi2 = psi_list1[ii], psi_list2[ii]
        vlsvs.append(VoilaLsv(bins, lsv_graphic=lsv_graphic, psi1=psi1, psi2=psi2))

    pickle.dump(VoilaInput(vlsvs, meta_info), open(pickle_path, 'w'))


def __load_execution_chunk(grp1_fnames, grp2_fnames, prior_path):


    matched_files = [None] * len(grp1_fnames)
    for idx, fname in enumerate(grp1_fnames):
        matched_files[idx] =  pickle.load(open(fname))
    filtered_lsv1 = majiq_filter.quantifiable_in_group(matched_files, 1, 1, None, None)

    matched_files = [None] * len(grp2_fnames)
    for idx, fname in enumerate(grp2_fnames):
        matched_files[idx] = pickle.load(open(fname))
    filtered_lsv2 = majiq_filter.quantifiable_in_group(matched_files, 1, 1, None, None)

    matched_lsv, matched_info = majiq_filter.lsv_intersection(filtered_lsv1, filtered_lsv2, bycol=True)

    print "After intersection:  %d/(%d, %d)" % (len(matched_info), len(filtered_lsv1[0]),
                                                      len(filtered_lsv2[0]))

    psi_space = np.linspace(0, 1 - BINSIZE, num=NUMBINS) + BINSIZE / 2
    prior = pickle.load(open(prior_path))

    return np.array(matched_lsv[0]), np.array(matched_lsv[1]), matched_info, psi_space, prior

def majiq_dpsi(fname1, fname2, delta_prior_path, num_exp, m_samples=100, debug=False):

    lsv_samples1, lsv_samples2, info, psi_space, prior_matrix = __load_execution_chunk( fname1, fname2,
                                                                                        delta_prior_path)

    nbins = len(psi_space)
    print "Calculating deltas..."

    post_matrix = []
    new_info = []
    ones_n = np.ones(shape=(1, nbins), dtype=np.float)

    # pickle.dump([lsv_samples1, info], open('./lsv_binomproblem.pkl', 'w+b'))
    posterior_psi1 = []
    posterior_psi2 = []
    #print lsv_samples1
    if debug:
        info = info[:100]
    for lidx, lsv_info in enumerate(info):
        num_ways = len(lsv_samples1[lidx][0])
        if lidx % 50 == 0:
            print "Event %d ..." % lidx,
            sys.stdout.flush()

        alpha_prior, beta_prior = majiq_psi.__get_prior_params(lsv_info, num_ways)
        if 'i' in lsv_info[2]:
            prior_idx = 1
        else:
            prior_idx = 0
        post_matrix.append([])
        new_info.append(lsv_info)

        posterior_psi1.append([])
        posterior_psi2.append([])
        #print lidx, lsv_samples1.shape
        psi1 = lsv_samples1[lidx, :]
        psi2 = lsv_samples2[lidx, :]

        for p_idx in xrange(num_ways):

            posterior = np.zeros(shape=(nbins, nbins), dtype=np.float)
            post_psi1 = np.zeros(shape=nbins, dtype=np.float)
            post_psi2 = np.zeros(shape=nbins, dtype=np.float)
            for m in xrange(m_samples):
                # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                junc = np.array([psi1[xx][p_idx][m] for xx in xrange(num_exp[0])])
                all_sample = np.array([psi1[xx][yy][m].sum() for xx in xrange(num_exp[0]) for yy in xrange(num_ways)])
                data_given_psi1 = np.log(majiq_psi.prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                    alpha_prior[p_idx], beta_prior[p_idx]))

                psi_v1 = data_given_psi1.reshape(nbins, -1)
                post_psi1 += (data_given_psi1 - scipy.misc.logsumexp(data_given_psi1))

                junc = np.array([psi2[xx][p_idx][m] for xx in xrange(num_exp[1])])
                all_sample = np.array([psi2[xx][yy][m].sum() for xx in xrange(num_exp[1]) for yy in xrange(num_ways)])

                data_given_psi2 = np.log(majiq_psi.prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                    alpha_prior[p_idx], beta_prior[p_idx]))
                post_psi2 += (data_given_psi2 - scipy.misc.logsumexp(data_given_psi2))
                psi_v2 = data_given_psi2.reshape(-1, nbins)

                A = (psi_v1 * ones_n + psi_v2 * ones_n.T) + np.log(prior_matrix[prior_idx])
                posterior += np.exp(A - scipy.misc.logsumexp(A))

            post_matrix[-1].append(posterior / m_samples)
            posterior_psi1[-1].append(post_psi1 / m_samples)
            posterior_psi2[-1].append(post_psi2 / m_samples)
            if num_ways == 2:
                break

    return post_matrix, new_info, posterior_psi1, posterior_psi2


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-grp1', dest="files1", nargs='+', required=True)
    parser.add_argument('-grp2', dest="files2", nargs='+', required=True)
    parser.add_argument('-prior', dest="prior_file", required=True)
    args = parser.parse_args()

    num_exp = (len(args.files1), len(args.files2))
    KK = majiq_dpsi(args.files1, args.files2, args.prior_file, num_exp, m_samples=100, debug=False)
    print "Wrapping up"
    _dump_lsvs_voila('majiq_vals.pickle', KK[0], KK[1], None, KK[2], KK[3])
