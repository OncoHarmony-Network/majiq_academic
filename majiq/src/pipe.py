import sys
import os
from multiprocessing import current_process

import scipy.misc
import numpy as np

import majiq.src.io_utils
from majiq.src.psi import prob_data_sample_given_psi, __get_prior_params
from majiq.src.utils.utils import get_logger
import majiq.src.io as majiq_io
import majiq.src.sample as majiq_sample


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

    sys.stdout.flush()
    thread_logger.info("[Th %s]: Saving DeltaPSI..." % chunk)
    majiq_io.dump_bin_file(results, "%s/%s_th%s.%s.pickle" % (tempdir, name, chunk, func.__name__))


def __load_execution_chunk(filename, delta=None):
    l_vals = majiq.src.io_utils.load_bin_file(filename)
    if not delta is None:
        prior = majiq.src.io_utils.load_bin_file(delta)
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

    #
    # #TODO: DELETE THIS ONLY TEMP
    # print "IF YOU READ THIS, your code needs to be changed, ask Jordi"
    #
    # num_exp = len(matched_lsv[0])
    # print num_exp
    # dir = './to_scott'
    # if not os.path.isdir(dir):
    #     os.mkdir(dir)
    # for ii in range(num_exp):
    #     res = []
    #     for lidx, lsv_info in enumerate(info):
    #         res.append(lsv_samples[lidx, ii])
    #
    #
    #     outfp = open('%s/%s.%d.pickle' % (dir, conf['name'], ii), 'w+bcott')
    #     pickle.dump([res, info], outfp)
    #     outfp.close()
    # exit()


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
                all_sample = np.array([psi[xx][yy][m].sum() for xx in xrange(num_exp) for yy in xrange(num_ways)])
                data_given_psi = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                   alpha_prior[p_idx], beta_prior[p_idx]))
                # normalizing
                posterior += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))

            post_psi[-1].append(posterior / conf['m'])
            if num_ways == 2:
                break


    return post_psi, new_info


# def deltapsi(matched_lsv, info, num_exp, conf, prior_matrix,  fitfunc, psi_space, logger):
