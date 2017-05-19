import h5py
import scipy.misc
import sys
import multiprocessing as mp
import traceback

from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
import majiq.src.utils as majiq_utils
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params, samples_from_psi, bootstrap_samples_calculation
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager

from voila.io_voila import Voila
from voila.vlsv import Het
from majiq.src.stats import operator, all_stats


def het_quantification(args_vals):
    list_of_lsv, chnk = args_vals
    logger = majiq_utils.get_logger("%s/%s.majiq.log" % (process_conf.outDir, chnk),
                                    silent=process_conf.silent, debug=process_conf.debug)
    try:
        logger.info("Quantifying LSVs PSI.. %s" % chnk)
        num_exp = [len(process_conf.files1), len(process_conf.files2)]
        files = [process_conf.files1, process_conf.files2]
        fitfunc = [None, None]
        lsvs = [None, None]

        for grp_idx in range(2):
            lsvs[grp_idx], fitfunc[grp_idx] = majiq_io.get_extract_lsv_list(list_of_lsv, files[grp_idx])

        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print "Event %d ..." % lidx
                sys.stdout.flush()

            lsv_het = Het()

            samps = [None, None]
            for grp_idx in range(2):
                quant_lsv = lsvs[grp_idx][lidx]
                lsv_type = quant_lsv.type
                boots = bootstrap_samples_calculation(quant_lsv, n_replica=num_exp[grp_idx],
                                                      name=process_conf.names[grp_idx],
                                                      outdir=process_conf.outDir,
                                                      nbins=process_conf.nbins,
                                                      store_bootsamples=False,
                                                      lock_array=process_conf.lock_per_file,
                                                      fitfunc_r=fitfunc[grp_idx],
                                                      m_samples=process_conf.m,
                                                      k_positions=process_conf.k,
                                                      discardzeros=process_conf.discardzeros,
                                                      trimborder=process_conf.trimborder,
                                                      debug=process_conf.debug)

                num_ways = boots.shape[1]
                alpha_prior, beta_prior = get_prior_params(lsv_type, num_ways)
                mu_psi = np.zeros(shape=(num_exp[grp_idx], num_ways))
                mean_psi = np.zeros(shape=(num_ways, process_conf.nbins), dtype=np.float)
                samps[grp_idx] = np.zeros(shape=(num_exp[grp_idx], num_ways, process_conf.nsamples))
                for exp in xrange(num_exp[grp_idx]):
                    all_sample = boots[exp].sum(axis=0)
                    for p_idx in xrange(num_ways):
                        alpha_0 = alpha_prior[p_idx]
                        beta_0 = beta_prior[p_idx]
                        post_psi = np.zeros(shape=process_conf.nbins, dtype=np.float)
                        for m in xrange(process_conf.m):
                            junc = boots[exp, p_idx, m]
                            data_given_psi = np.log(prob_data_sample_given_psi(junc, all_sample[m],
                                                                               process_conf.nbins,
                                                                               alpha_0, beta_0))
                            post_psi += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))
                            mu_psi[exp, p_idx] += float(junc + alpha_0) / (all_sample[m] + alpha_0 + beta_0)

                        post_psi /= process_conf.m
                        mean_psi[p_idx] += post_psi
                        mu_psi[exp, p_idx] /= process_conf.m

                        if process_conf.nsamples == 1:
                            samps[grp_idx][exp, p_idx, 0] = mu_psi
                        else:
                            samps[grp_idx][exp, p_idx, :] = samples_from_psi(post_psi, mu_psi[exp, p_idx], process_conf.vwindow,
                                                                             process_conf.nsamples,
                                                                             process_conf.nbins)

                mean_psi /= num_exp[grp_idx]
                lsv_het.add_group(mu_psi, mean_psi)
                print grp_idx, mu_psi

            out_stats = do_test_stats(samps, process_conf.stats, process_conf.minsamps)
            for stat_idx in xrange(out_stats.shape[1]):
                lsv_het.add_junction_stats(out_stats[:, stat_idx])

                # if num_ways == 2:
                #     break

            qm = QueueMessage(QUEUE_MESSAGE_HETER_DELTAPSI, (lsv_het, lsv_id), chnk)
            process_conf.queue.put(qm, block=True)

        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        process_conf.queue.put(qm, block=True)
        process_conf.lock[chnk].acquire()
        process_conf.lock[chnk].release()

    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        logger.error(e)
        raise()


def do_test_stats(insamps, stats, minsamps):
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

    cts = [insamps[0].shape[0], insamps[1].shape[0]]
    # lsv_id, lsv_type = lsv_attrs
    # samps = []
    # cts = []
    # for grp in samples.values():
    #     cts.append(len(grp))
    #     samps.append(grp)
    samps = np.concatenate(insamps, axis=0)
    labels = np.concatenate((np.zeros(cts[0]), np.ones(cts[1])), axis=0).astype(int)
    nfiles, njuncs, nsamples = samps.shape
    outstats = np.zeros(shape=(njuncs, len(stats)))
    if all([nsamps > minsamps for nsamps in cts]):
        for jn_idx in range(njuncs):
            cur_output = np.zeros(shape=(len(stats), nsamples))
            for samp_idx in range(nsamples):
                csamps = samps[:, jn_idx, samp_idx]
                clabels = labels[csamps != -1]
                npos = clabels.sum()
                nneg = clabels.size - npos
                if npos < minsamps or nneg < minsamps:
                    continue
                csamps = csamps[csamps != -1]
                asort = csamps.argsort()
                csamps = csamps[asort]
                clabels = clabels[asort]
                for stat_idx, stat_name in enumerate(stats):
                    cur_output[stat_idx, samp_idx] = operator[stat_name].operator(csamps, clabels)

            outstats[jn_idx, :] = np.percentile(cur_output, 95, axis=1)

    return outstats


def calc_independent(args):
    pipeline_run(independent(args))


class independent(BasicPipeline):
    def run(self):
        self.independent()

    def independent(self):
        majiq_utils.create_if_not_exists(self.logger_path)
        logger = majiq_utils.get_logger("%s/deltapsi_majiq.log" % self.logger_path, silent=self.silent,
                                        debug=self.debug)

        logger.info("Majiq deltapsi v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("GROUP1: %s" % self.files1)
        logger.info("GROUP2: %s" % self.files2)
        self.nbins = 40

        try:
            for stats_name in self.stats:
                module_ = __import__('majiq.src.stats.'+stats_name.lower(), fromlist=stats_name.title())
                class_ = getattr(module_, stats_name.title())
                operator[stats_name] = class_()
        except ImportError as i_err:
            logger.error("The %s statistic is not one of the available statistics, "
                         "in  [ %s ]" % (stats_name, ' | '.join(all_stats)))
            return

        lsv_dict1, lsv_types1, lsv_summarized1, meta1, lsv_dict_graph1 = majiq_io.extract_lsv_summary(self.files1)
        list_of_lsv_group1 = majiq_filter.merge_files_hdf5(lsv_dict1, lsv_summarized1, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=logger)

        lsv_dict2, lsv_types2, lsv_summarized2, meta2, lsv_dict_graph2 = majiq_io.extract_lsv_summary(self.files2)
        list_of_lsv_group2 = majiq_filter.merge_files_hdf5(lsv_dict2, lsv_summarized2, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=logger)

        list_of_lsv = list(set(list_of_lsv_group1).intersection(set(list_of_lsv_group2)))
        lchnksize = max(len(list_of_lsv) / self.nthreads, 1) + 1
        if len(list_of_lsv) > 0:
            q = mp.Queue()
            lock_arr = [mp.Lock() for xx in range(self.nthreads)]
            [xx.acquire() for xx in lock_arr]
            pool = mp.Pool(processes=self.nthreads, initializer=process_conf,
                           initargs=[self, q, lock_arr, None],
                           maxtasksperchild=1)

            pool.map_async(het_quantification, majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            with Voila(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
                out_h5p.add_genome(meta1['genome'])
                out_h5p.add_experiments(self.names[0], experiment_names=meta1['experiments'])
                out_h5p.add_experiments(self.names[1], experiment_names=meta2['experiments'])
                out_h5p.add_stat_names(self.stats)

                in_h5p = h5py.File(self.files1[0], 'r')
                queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=self.nthreads, logger=logger,
                              list_of_lsv_graphics=lsv_dict_graph1)
                in_h5p.close()
            pool.join()
