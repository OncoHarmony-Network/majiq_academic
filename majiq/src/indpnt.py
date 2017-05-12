import h5py
import numpy as np
import scipy.misc
import sys
import multiprocessing as mp
import traceback

from majiq.src.basic_pipeline import BasicPipeline, pipeline_run, bootstrap_samples_with_divs
import majiq.src.utils as majiq_utils
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
from majiq.src.io_utils import dump_bin_file, load_bin_file
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params, samples_from_psi
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, quantification_init, queue_manager


def het_quantification(args_vals):
    list_of_lsv, chnk = args_vals
    logger = majiq_utils.get_logger("%s/%s.majiq.log" % (quantification_init.output, chnk),
                                        silent=quantification_init.silent, debug=quantification_init.debug)
    try:
        logger.info("Quantifying LSVs PSI.. %s" % chnk)
        num_exp = [len(quantification_init.files[0]), len(quantification_init.files[1])]

        fitfunc = [None, None]
        lsvs = [None, None]
        for grp_idx in range(2):
            lsvs[grp_idx], fitfunc[grp_idx] = majiq_io.get_extract_lsv_list(list_of_lsv,
                                                                            quantification_init.files[grp_idx])

        prior_matrix = np.array(load_bin_file(get_prior_matrix_filename(quantification_init.output,
                                                                        quantification_init.names)))

        ones_n = np.ones(shape=(1, quantification_init.nbins), dtype=np.float)

        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print "Event %d ..." % lidx
                sys.stdout.flush()
            lsv_samples = [None, None]
            prior_idx = 1 if 'i' in lsv_type else 0
            num_ways = psi1.shape[1]
            alpha_prior, beta_prior = get_prior_params(lsv_type, num_ways)
            for grp_idx in range(2):

                quant_lsv = lsvs[grp_idx][lidx]
                lsv_type = quant_lsv.type
                boots = bootstrap_samples_calculation(quant_lsv, n_replica=num_exp[grp_idx],
                                                                     name=quantification_init.names[grp_idx],
                                                                     outdir=quantification_init.output,
                                                                     nbins=quantification_init.nbins,
                                                                     store_bootsamples=quantification_init.boots,
                                                                     lock_array=quantification_init.lock_per_file,
                                                                     fitfunc_r=fitfunc[grp_idx],
                                                                     m_samples=quantification_init.m,
                                                                     k_positions=quantification_init.k,
                                                                     discardzeros=quantification_init.discardzeros,
                                                                     trimborder=quantification_init.trimborder,
                                                                     debug=quantification_init.debug)

                all_sample = boots.sum(axis=0)
                post_psi = np.zeros(shape=quantification_init.nbins, dtype=np.float)
                for exp in xrange(num_exp[grp_idx]):
                    for p_idx in xrange(num_ways):
                        alpha_0 = alpha_prior[p_idx]
                        beta_0 = beta_prior[p_idx]
                        mu_psi = 0
                        for m in xrange(quantification_init.m):
                            junc = boots[p_idx, m]
                            data_given_psi = np.log(prob_data_sample_given_psi(junc, all_sample[m],
                                                                                quantification_init.nbins,
                                                                                alpha_0, beta_0))
                            post_psi += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))
                            mu_psi += float(junc + alpha_0) / (all_sample[m] + alpha_0 + beta_0)

                        post_psi /= quantification_init.m
                        mu_psi /= quantification_init.m

                        if quantification_init.nsamples == 1:
                            samps[grp_idx][p_idx, 0] = mu_psi
                        else:
                            samps[grp_idx][p_idx, :] = samples_from_psi(post_psi, mu_psi, quantification_init.vwindow,
                                                               quantification_init.nsamples, quantification_init.nbins)

                do_test_stats(samps, stats, minsamps)

                if num_ways == 2:
                    break

            qm = QueueMessage(QUEUE_MESSAGE_DELTAPSI_RESULT, (post_matrix, posterior_psi1, posterior_psi2,
                                                              mu_psi1, mu_psi2, lsv_id), chnk)
            quantification_init.queue.put(qm, block=True)

        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        quantification_init.queue.put(qm, block=True)
        quantification_init.lock[chnk].acquire()
        quantification_init.lock[chnk].release()

    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


def do_test_stats(samps, stats, minsamps):
#def do_test_stats(lsv_attrs, samples, stats, minsamps):
    """
    Log P values for the given set of statistical tests
    :param lsv_attrs: 2-tuple (lsv_id, lsv_type)
    :param samples: 3d numpy array with dimensions (nfiles, njunc, nsamps)
    :param stats: iterable of test statistics to perform
    :param minsamps: minimum number of samples to perform the test
    :return: 3d numpy array with dimensions
             (njunc, nsamps, len(stats)) containing log P values
    """

    cts = [samps[0].shape[0], samps[1].shape[0]]
    # lsv_id, lsv_type = lsv_attrs
    # samps = []
    # cts = []
    # for grp in samples.values():
    #     cts.append(len(grp))
    #     samps.append(grp)
    samps = np.vstack(samps)
    labels = np.hstack((np.zeros(cts[0]), np.ones(cts[1]))).astype(int)
    nfiles, njuncs, nsamples = samps.shape
    outstats = np.zeros((njuncs, nsamples, len(stats)))
    if all([nsamps > minsamps for nsamps in cts]):
        for jn_idx in range(njuncs):
            for samp_idx in range(nsamples):
                csamps = samps[:, jn_idx, samp_idx]
                clabels = labels[csamps != -1]
                npos = clabels.sum()
                nneg = clabels.size - npos
                if npos < minsamps or nneg < minsamps:
                    continue
                csamps = csamps[csamps != -1]
                nsamps = csamps[clabels == 0]
                psamps = csamps[clabels == 1]
                asort = csamps.argsort()
                ssamps = csamps[asort]
                slabels = labels[asort]
                for stat_idx, stat_name in enumerate(stats):
                    stat_name = stat_name.title()
                    #op, flg = operators[stat_name]
                    if flg:
                        cur_output = op(ssamps, slabels)
                    else:
                        cur_output = op(nsamps, psamps)
                    outstats[jn_idx, samp_idx, stat_idx] = cur_output
    #return lsv_id, lsv_type, outstats
    return outstats


def independent(args):
    pipeline_run(independent(args))


class independent(BasicPipeline):
    def run(self):
        self.independent()

    def independent(self):
        majiq_utils.create_if_not_exists(self.logger_path)
        self.logger = majiq_utils.get_logger("%s/deltapsi_majiq.log" % self.logger_path, silent=self.silent,
                                             debug=self.debug)

        self.logger.info("Majiq deltapsi v%s" % VERSION)
        self.logger.info("Command: %s" % " ".join(sys.argv))
        self.logger.info("GROUP1: %s" % self.files1)
        self.logger.info("GROUP2: %s" % self.files2)
        self.nbins = 20

        files = [self.files1, self.files2]
        file_locks = None
        if self.export_boots:
            file_locks = [[mp.Lock() for xx in self.files1], [mp.Lock() for xx in self.files2]]

        # lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()
        lock_arr = [mp.Lock() for xx in range(self.nthreads)]

        lsv_dict1, lsv_types1, lsv_summarized1, meta1 = majiq_io.extract_lsv_summary(self.files1)
        list_of_lsv_group1 = majiq_filter.merge_files_hdf5(lsv_dict1, lsv_summarized1, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=self.logger)

        lsv_dict2, lsv_types2, lsv_summarized2, meta2 = majiq_io.extract_lsv_summary(self.files2)
        list_of_lsv_group2 = majiq_filter.merge_files_hdf5(lsv_dict2, lsv_summarized2, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=self.logger)

        list_of_lsv = list(set(list_of_lsv_group1).intersection(set(list_of_lsv_group2)))

        if len(list_of_lsv) > 0:
            [xx.acquire() for xx in lock_arr]
            pool = mp.Pool(processes=self.nthreads, initializer=quantification_init,
                           initargs=[q, lock_arr, self.output, self.name, self.silent, self.debug, self.nbins, self.m,
                                     self.k, self.discardzeros, self.trimborder, self.files, self.only_boots, weights,
                                     file_locks],
                           maxtasksperchild=1)

            pool.map_async(het_quantification, majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            with Voila(get_quantifier_voila_filename(self.output, self.name), 'w') as out_h5p:
                out_h5p.add_metainfo(meta['genome'], self.name, meta['experiments'])
                in_h5p = h5py.File(self.files[0], 'r')
                queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=self.nthreads, logger=self.logger)
                in_h5p.close()

            pool.join()
