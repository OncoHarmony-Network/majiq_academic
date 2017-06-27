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
import majiq.src.sample as majiq_sample
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params, gen_prior_matrix, bootstrap_samples_calculation
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager

from voila.io_voila import Voila
from voila.constants import ANALYSIS_DELTAPSI
import collections


def deltapsi(args):
    return pipeline_run(DeltaPsi(args))


def deltapsi_quantification(args_vals):
    try:
        list_of_lsv, chnk = args_vals
        logger = majiq_utils.get_logger("%s/%s.majiq.log" % (process_conf.outDir, chnk),
                                        silent=process_conf.silent, debug=process_conf.debug)

        logger.info("Quantifying LSVs PSI.. %s" % chnk)
        num_exp = [len(process_conf.files1), len(process_conf.files1)]

        f_list = [None, None]

        f_list[0] = majiq_io.get_extract_lsv_list(list_of_lsv, process_conf.files1)
        f_list[1] = majiq_io.get_extract_lsv_list(list_of_lsv, process_conf.files2)

        prior_matrix = np.array(load_bin_file(get_prior_matrix_filename(process_conf.outDir,
                                                                        process_conf.names)))

        ones_n = np.ones(shape=(1, process_conf.nbins), dtype=np.float)

        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print("Event %d ..." % lidx)
                sys.stdout.flush()
            lsv_samples = [None, None]

            for grp_idx in range(2):
                lsv_samples[grp_idx] = f_list[grp_idx][lidx].coverage
                lsv_type = f_list[grp_idx][lidx].type

            prior_idx = 1 if 'i' in lsv_type else 0

            psi1, psi2 = [np.array(xx) for xx in lsv_samples]
            msamples = psi1.shape[2]
            del lsv_samples

            post_matrix = []
            posterior_psi1 = []
            posterior_psi2 = []
            num_ways = psi1.shape[1]
            mu_psi1 = []
            mu_psi2 = []
            alpha_prior, beta_prior = get_prior_params(lsv_type, num_ways)
            for p_idx in range(num_ways):

                posterior = np.zeros(shape=(process_conf.nbins, process_conf.nbins), dtype=np.float)
                post_psi1 = np.zeros(shape=process_conf.nbins, dtype=np.float)
                post_psi2 = np.zeros(shape=process_conf.nbins, dtype=np.float)
                mu_psi1_m = []
                mu_psi2_m = []
                alpha_0 = alpha_prior[p_idx]
                beta_0 = beta_prior[p_idx]
                for m in range(msamples):
                    # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                    # junc = [psi1[xx][p_idx][m] for xx in xrange(num_exp[0])]
                    # junc = np.array(junc)
                    junc = psi1[:, p_idx, m]
                    all_sample = [psi1[xx][yy][m].sum() for xx in range(num_exp[0]) for yy in range(num_ways)]
                    all_sample = np.array(all_sample)
                    mu_psi1_m.append(float(junc.sum() + alpha_0) / (all_sample.sum() + alpha_0 + beta_0))

                    data_given_psi1 = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(),
                                                                        process_conf.nbins,
                                                                        alpha_0, beta_0))

                    psi_v1 = data_given_psi1.reshape(process_conf.nbins, -1)
                    post_psi1 += np.exp(data_given_psi1 - scipy.misc.logsumexp(data_given_psi1))

                    # junc = [psi2[xx][p_idx][m] for xx in xrange(num_exp[1])]
                    # junc = np.array(junc)
                    junc = psi2[:, p_idx, m]
                    all_sample = [psi2[xx][yy][m].sum() for xx in range(num_exp[1]) for yy in range(num_ways)]
                    all_sample = np.array(all_sample)
                    mu_psi2_m.append(float(junc.sum() + alpha_0) / (all_sample.sum() + alpha_0 + beta_0))
                    data_given_psi2 = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(),
                                                                        process_conf.nbins,
                                                                        alpha_0, beta_0))
                    post_psi2 += np.exp(data_given_psi2 - scipy.misc.logsumexp(data_given_psi2))
                    psi_v2 = data_given_psi2.reshape(-1, process_conf.nbins)

                    A = (psi_v1 * ones_n + psi_v2 * ones_n.T) + np.log(prior_matrix[prior_idx])

                    posterior += np.exp(A - scipy.misc.logsumexp(A))
                mu_psi1.append(np.median(mu_psi1_m))
                mu_psi2.append(np.median(mu_psi2_m))
                post_matrix.append(posterior / msamples)
                posterior_psi1.append(post_psi1 / msamples)
                posterior_psi2.append(post_psi2 / msamples)
                if num_ways == 2:
                    break

            qm = QueueMessage(QUEUE_MESSAGE_DELTAPSI_RESULT, (post_matrix, posterior_psi1, posterior_psi2,
                                                              mu_psi1, mu_psi2, lsv_id), chnk)
            process_conf.queue.put(qm, block=True)

        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        process_conf.queue.put(qm, block=True)
        process_conf.lock[chnk].acquire()
        process_conf.lock[chnk].release()

    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


conf = collections.namedtuple('conf', 'name output_dir silent debug markstacks')
prior_conf = collections.namedtuple('conf', 'iter plotpath breakiter names binsize')


class DeltaPsi(BasicPipeline):

    def run(self):
        self.deltapsi()

    def deltapsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """

        majiq_utils.create_if_not_exists(self.logger_path)
        self.logger = majiq_utils.get_logger("%s/deltapsi_majiq.log" % self.logger_path, silent=self.silent,
                                             debug=self.debug)

        self.logger.info("Majiq deltapsi v%s" % VERSION)
        self.logger.info("Command: %s" % " ".join(sys.argv))
        self.logger.info("GROUP1: %s" % self.files1)
        self.logger.info("GROUP2: %s" % self.files2)
        self.nbins = 20

        files = [self.files1, self.files2]

        # lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()
        lock_arr = [mp.Lock() for xx in range(self.nthreads)]

        lsv_dict1, lsv_types1, lsv_summarized1, meta1, lsv_dict_graph1 = majiq_io.extract_lsv_summary(self.files1)
        list_of_lsv_group1 = majiq_filter.merge_files_hdf5(lsv_dict1, lsv_summarized1, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=self.logger)

        print (len(list_of_lsv_group1))

        lsv_dict2, lsv_types2, lsv_summarized2, meta2, lsv_dict_graph2 = majiq_io.extract_lsv_summary(self.files2)
        list_of_lsv_group2 = majiq_filter.merge_files_hdf5(lsv_dict2, lsv_summarized2, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=self.logger)
        print (len(list_of_lsv_group2))
        list_of_lsv = list(set(list_of_lsv_group1).intersection(set(list_of_lsv_group2)))


        psi_space, prior_matrix = gen_prior_matrix(lsv_dict1, lsv_summarized1, lsv_dict2, lsv_summarized2, lsv_types1,
                                                   self.outDir, prior_conf(self.iter, self.plotpath, self.breakiter,
                                                                           self.names, self.binsize),
                                                   numbins=self.nbins, defaultprior=self.default_prior,
                                                   minpercent=self.min_exp, logger=self.logger)

        self.logger.info("Saving prior matrix for %s..." % self.names)
        dump_bin_file(prior_matrix, get_prior_matrix_filename(self.outDir, self.names))

        lchnksize = max(len(list_of_lsv)/self.nthreads, 1) + 1
        print (len(list_of_lsv))
        weights = [None, None]
        for grp_idx, fil in enumerate(files):
            #TODO: FIX WEIGHTS, maybe just return rho
            weights[grp_idx] = self.calc_weights(self.weights[grp_idx], fil, list_of_lsv, lock_arr, lchnksize, q,
                                                 self.names[grp_idx])

        pool = mp.Pool(processes=self.nthreads, initializer=process_conf, initargs=[self, q, lock_arr, weights],
                       maxtasksperchild=1)
        [xx.acquire() for xx in lock_arr]

        if len(list_of_lsv) > 0:
            pool.map_async(deltapsi_quantification,
                           majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            with Voila(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
                #out_h5p.add_genome(meta1['genome'])
                out_h5p.set_analysis_type(ANALYSIS_DELTAPSI)
                out_h5p.add_experiments(group_name=self.names[0], experiment_names=meta1['experiments'])
                out_h5p.add_experiments(group_name=self.names[1], experiment_names=meta2['experiments'])
                # out_h5p.add_metainfo(meta1['genome'], group1=self.names[0], experiments1=meta1['experiments'],
                #                      group2=self.names[1], experiments2=meta2['experiments'])

                in_h5p = h5py.File(files[0][0], 'r')
                queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=self.nthreads, logger=self.logger,
                              list_of_lsv_graphics=lsv_dict_graph1)
                in_h5p.close()
            pool.join()

        self.logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                          self.names[1],
                                                                                                          self.outDir))
        self.logger.info("Alakazam! Done.")
