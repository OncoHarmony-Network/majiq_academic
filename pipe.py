import analysis.psi as majiq_psi
import scipy.misc
import random
import numpy as np
from scipy.stats import pearsonr, binom_test
from scipy.stats import binom, beta
from grimoire.utils.utils import create_if_not_exists, get_logger
import sys
import os
from multiprocessing import Pool, Manager, current_process
import analysis.sample as majiq_sample 
import pickle


def parallel_lsv_child_calculation(func, args, info, tempdir, name, chunk):

    try:
        if not os.path.isdir(tempdir):
                os.mkdir(tempdir)
        thread_logger = get_logger("%s/majiq.w%s.log" % (tempdir, chunk), silent=False)
        thread_logger.info("[Th %s]: START child,%s" % (chunk, current_process().name))
        thread_logger.info('[Th %s]: Filtering ...' % chunk)

        args.append(thread_logger)
#        post_matrix, new_info = model2( *args )
        post_matrix, new_info = func(*args)

        print "%s/%s_th%s.deltapsi.pickle" % (tempdir, name, chunk)
        sys.stdout.flush()
        thread_logger.info("[Th %s]: Saving DeltaPSI..." % chunk)
        output = open("%s/%s_th%s.%s.pickle" % (tempdir, name, chunk, func.__name__), 'w')
        pickle.dump([post_matrix, new_info], output)

    except Exception as e:
        print "%s" % sys.exc_traceback.tb_lineno, e
        sys.stdout.flush()

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


def prob_data_sample_given_psi(sample, all_sample, psi_space, nbins):
    bsize = bsize = 1.0/float(nbins)
    if all_sample <= 4000:
        #bin_test = [binom_test(sample, all_sample, p=x) for x in psi_space]
        bin_test = [binom.pmf(sample, all_sample, p=x) for x in psi_space]
    else:
        psi_border = np.arange(0, 1.01, bsize)
        notsample = all_sample - sample
        bincdf = [beta.cdf(xx, sample + 0.5, notsample + 0.5) for xx in psi_border]

        bin_test = []
        for x in xrange(nbins):
            val = bincdf[x] - bincdf[x+1]
            bin_test.append(val)

    bin_test = np.array(bin_test) + 1e-10
    return bin_test


def calcpsi(matched_lsv, info, num_exp, conf, fitfunc, logger):

    try:

        bsize = 1.0/float(conf['nbins'])
        psi_space = np.arange(0 + bsize/2, 1, bsize)
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
                                                                           debug=conf['debug'],
                                                                           Nz=conf['nz'])

                lsv_samples[lidx, eidx] = s_lsv

        nbins = conf['nbins']
        logger.info("Calculating psis...")

        post_psi = []
        new_info = []
        for lidx, lsv_info in enumerate(info):
            num_ways = float(len(lsv_samples[lidx, 0]))
            if lidx % 50 == 0:
                print "Event %d ..." % lidx,
                sys.stdout.flush()

            alpha = 1.0/num_ways
            beta = (num_ways-1.0) / num_ways
            prior = majiq_psi.dirichlet_pdf(np.array([psi_space, 1-psi_space]).T, np.array([alpha, beta]))

            post_psi.append([])
            new_info.append(lsv_info)
            for p_idx in xrange(int(num_ways)):
                posterior = np.zeros(shape=nbins, dtype=np.float)
                psi = np.zeros(shape=nbins, dtype=np.float)
                for m in xrange(conf['m']):
                    # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                    data_given_psi = np.zeros(shape=(num_exp, nbins), dtype = np.float)
                    for exp_idx in xrange(num_exp):
                        junc = lsv_samples[lidx, exp_idx][p_idx][m]
                        all_psi = np.array([xx[m] for xx in lsv_samples[lidx, exp_idx]])
                        data_given_psi[exp_idx] = np.log(prob_data_sample_given_psi(junc, all_psi.sum(),
                                                                                    psi_space, conf['nbins']))
                    psi += data_given_psi.sum(axis=0) + np.log(prior)
                    # normalizing
                    posterior += np.exp(psi - scipy.misc.logsumexp(psi))

                post_psi[-1].append(posterior / conf['m'])
                if num_ways == 2:
                    break

    except Exception as e:
        post_psi = []
        new_info = []
        print "%s" % sys.exc_traceback.tb_lineno, e
        sys.stdout.flush()

    return post_psi, new_info


def model2(matched_lsv, info, num_exp, conf, prior_matrix,  fitfunc, psi_space, logger):

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
                                                                           debug=conf['debug'],
                                                                           Nz=conf['nz'])
                if grp_idx == 0:
                    lsv_samples1[lidx, eidx] = s_lsv
                else:
                    lsv_samples2[lidx, eidx] = s_lsv

    nbins = len(psi_space)
    logger.info("Calculating deltas...")

    post_matrix = []
    new_info = []
    ones_n = np.ones(shape=(1, nbins), dtype = np.float)
    for lidx, lsv_info in enumerate(info):
        num_ways = len(lsv_samples1[lidx][0])
        if lidx % 50 == 0:
            print "Event %d ..." % lidx,
            sys.stdout.flush()

        post_matrix.append([])
        new_info.append(lsv_info)
        for p_idx in xrange(num_ways):

            posterior = np.zeros(shape=(nbins, nbins), dtype = np.float)
            for m in xrange(conf['m']):
                # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                data_given_psi1 = np.zeros(shape=(num_exp[0], nbins), dtype = np.float)
                for exp_idx in xrange(num_exp[0]):
                    psi1 = lsv_samples1[lidx, exp_idx][p_idx][m]
                    all_psi = np.array([xx[m] for xx in lsv_samples1[lidx, exp_idx]])
                    data_given_psi1[exp_idx] = np.log(prob_data_sample_given_psi(psi1, all_psi.sum(),
                                                                                 psi_space, conf['nbins']))
                V1 = data_given_psi1.sum(axis=0)
                V1 = V1.reshape(nbins, -1)

                data_given_psi2 = np.zeros(shape=(num_exp[1], nbins), dtype = np.float)
                for exp_idx in xrange(num_exp[1]):
                    psi2 = lsv_samples2[lidx, exp_idx][p_idx][m]
                    all_psi = np.array([xx[m] for xx in lsv_samples2[lidx, exp_idx]])
                    data_given_psi2[exp_idx] = np.log(prob_data_sample_given_psi(psi2, all_psi.sum(),
                                                                                 psi_space, conf['nbins']))
                V2 = data_given_psi2.sum(axis=0)
                V2 = V2.reshape(-1, nbins)

                A = (V1 * ones_n + V2 * ones_n.T) + np.log(prior_matrix)
                posterior += np.exp(A - scipy.misc.logsumexp(A))
            post_matrix[-1].append(posterior / conf['m'])
            if num_ways == 2:
                break

    return post_matrix, new_info
