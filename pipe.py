import random
import sys
import os
from multiprocessing import current_process
import pickle

import scipy.misc

import numpy as np
from scipy.stats import beta

from src.utils.utils import get_logger
import src.sample as majiq_sample


def parallel_lsv_child_calculation(func, args, tempdir, name, chunk):
    # try:
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)
    thread_logger = get_logger("%s/majiq.w%s.log" % (tempdir, chunk), silent=False)
    thread_logger.info("[Th %s]: START child,%s" % (chunk, current_process().name))
    thread_logger.info('[Th %s]: Filtering ...' % chunk)

    args.append(thread_logger)
    # post_matrix, new_info = model2( *args )
    results = func(*args)

    print "%s/%s_th%s.deltapsi.pickle" % (tempdir, name, chunk)
    sys.stdout.flush()
    thread_logger.info("[Th %s]: Saving DeltaPSI..." % chunk)
    output = open("%s/%s_th%s.%s.pickle" % (tempdir, name, chunk, func.__name__), 'w')
    pickle.dump(results, output)
    #
    # except Exception as e:
    #     print "%s" % sys.exc_traceback.tb_lineno, e
    #     sys.stdout.flush()

    return


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


def prob_data_sample_given_psi(sample, all_sample, nbins, alpha_prior, beta_prior):
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


def __load_execution_chunk(filename, delta=None):
    l_vals = pickle.load(open(filename))
    if not delta is None:
        prior = pickle.load(open(delta))
        l_vals.append(prior)

    return l_vals


def calcpsi(fname, conf, logger):
    # try:
    matched_lsv, info, num_exp, fitfunc = __load_execution_chunk(fname)
    #The center of the previous BINS. This is used to calculate the mean value of each bin.
    lsv_samples = np.zeros(shape=(len(info), num_exp), dtype=np.dtype('object'))
    logger.info("Bootstrapping for all samples...")
    for lidx, lsv_all in enumerate(matched_lsv):
        for eidx, lsv in enumerate(lsv_all):
            m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(junction_list=lsv,
                                                                       m=conf['m'],
                                                                       k=conf['k'],
                                                                       discardzeros=conf['discardzeros'],
                                                                       trimborder=conf['trimborder'],
                                                                       fitted_one_over_r=fitfunc[eidx],
                                                                       debug=conf['debug'])

            lsv_samples[lidx, eidx] = s_lsv

    nbins = conf['nbins']
    logger.info("Calculating psis...")

    post_psi = []
    new_info = []



    #TODO: DELETE THIS ONLY TEMP
    print "IF YOU READ THIS, your code needs to be changed, ask Jordi"

    num_exp = len(matched_lsv[0])
    print num_exp
    dir = './to_scott'
    if not os.path.isdir(dir):
        os.mkdir(dir)
    for ii in range(num_exp):
        res = []
        for lidx, lsv_info in enumerate(info):
            res.append(lsv_samples[lidx, ii])


        outfp = open('./%s/%s.%d.pickle' % (dir, conf['name'], ii))
        pickle.dumps([res, info], outfp)
        outfp.close()
    exit()



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

            for m in xrange(conf['m']):
                # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                junc = [psi[xx][p_idx][m] for xx in xrange(num_exp)]
                junc = np.array(junc)
                all_sample = [psi[xx][yy][m].sum() for xx in xrange(num_exp) for yy in xrange(num_ways)]
                all_sample = np.array(all_sample)
                data_given_psi = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                   alpha_prior[p_idx], beta_prior[p_idx]))

                # normalizing
                posterior += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))

            post_psi[-1].append(posterior / conf['m'])
            if num_ways == 2:
                break

                # except Exception as e:
                # post_psi = []
                #     new_info = []
                #     print "%s" % sys.exc_traceback.tb_lineno, e
                #     sys.stdout.flush()

    return post_psi, new_info


# def deltapsi(matched_lsv, info, num_exp, conf, prior_matrix,  fitfunc, psi_space, logger):
def deltapsi(fname, delta_prior_path, logger):
    matched_lsv, info, num_exp, conf, fitfunc, psi_space, prior_matrix = __load_execution_chunk(fname,
                                                                                                delta=delta_prior_path)

    lsv_samples1 = np.zeros(shape=(len(info), num_exp[0]), dtype=np.dtype('object'))
    lsv_samples2 = np.zeros(shape=(len(info), num_exp[1]), dtype=np.dtype('object'))

    logger.info("Bootstrapping for all samples...")
    for grp_idx, group in enumerate(matched_lsv):
        for lidx, lsv_all in enumerate(group):
            for eidx, lsv in enumerate(lsv_all):
                m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(junction_list=lsv,
                                                                           m=conf['m'],
                                                                           k=conf['k'],
                                                                           discardzeros=conf['discardzeros'],
                                                                           trimborder=conf['trimborder'],
                                                                           fitted_one_over_r=fitfunc[grp_idx][eidx],
                                                                           debug=conf['debug'])
                if grp_idx == 0:
                    lsv_samples1[lidx, eidx] = s_lsv
                else:
                    lsv_samples2[lidx, eidx] = s_lsv

    nbins = len(psi_space)
    logger.info("Calculating deltas...")

    post_matrix = []
    new_info = []
    ones_n = np.ones(shape=(1, nbins), dtype=np.float)

    # pickle.dump([lsv_samples1, info], open('./lsv_binomproblem.pkl', 'w+b'))
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
            for m in xrange(conf['m']):
                # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                junc = [psi1[xx][p_idx][m] for xx in xrange(num_exp[0])]
                junc = np.array(junc)
                all_sample = [psi1[xx][yy][m].sum() for xx in xrange(num_exp[0]) for yy in xrange(num_ways)]
                all_sample = np.array(all_sample)
                data_given_psi1 = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                    alpha_prior[p_idx], beta_prior[p_idx]))

                psi_v1 = data_given_psi1.reshape(nbins, -1)
                post_psi1 += np.exp(data_given_psi1 - scipy.misc.logsumexp(data_given_psi1))

                junc = [psi2[xx][p_idx][m] for xx in xrange(num_exp[1])]
                junc = np.array(junc)
                all_sample = [psi2[xx][yy][m].sum() for xx in xrange(num_exp[1]) for yy in xrange(num_ways)]
                all_sample = np.array(all_sample)
                data_given_psi2 = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                    alpha_prior[p_idx], beta_prior[p_idx]))
                post_psi2 += np.exp(data_given_psi2 - scipy.misc.logsumexp(data_given_psi2))
                psi_v2 = data_given_psi2.reshape(-1, nbins)

                A = (psi_v1 * ones_n + psi_v2 * ones_n.T) + np.log(prior_matrix[prior_idx])
                posterior += np.exp(A - scipy.misc.logsumexp(A))

            post_matrix[-1].append(posterior / conf['m'])
            posterior_psi1[-1].append(post_psi1 / conf['m'])
            posterior_psi2[-1].append(post_psi2 / conf['m'])
            if num_ways == 2:
                break

    return post_matrix, new_info, posterior_psi1, posterior_psi2
