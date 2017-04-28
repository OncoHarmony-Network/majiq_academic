import multiprocessing as mp
import sys
import traceback
import h5py
import numpy as np
import scipy.misc
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.sample as majiq_sample
import majiq.src.utils as majiq_utils
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, quantification_init, queue_manager
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params, bootstrap_samples_calculation
from voila.io_voila import Voila

import collections


################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    return pipeline_run(CalcPsi(args))


def psi_quantification(args_vals):
    try:
        #print len(args_vals), args_vals
        list_of_lsv, chnk = args_vals
        logger = majiq_utils.get_logger("%s/%s.majiq.log" % (quantification_init.output, chnk),
                                        silent=quantification_init.silent, debug=quantification_init.debug)

        logger.info("Quantifying LSVs PSI.. %s" % chnk)
        if quantification_init.boots:
            logger.info(".. Only Bootstrap enabled .. %s" % chnk)

        nbins = quantification_init.nbins
        num_exp = len(quantification_init.files)

        if not quantification_init.weights:
            f_list = majiq_io.open_bootstrap_samples(num_exp=num_exp, directory=quantification_init.output,
                                                     name=quantification_init.names)
        else:
            lsvs, fitfunc = majiq_io.get_extract_lsv_list(list_of_lsv, quantification_init.files)

        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print "Event %d ..." % lidx
                sys.stdout.flush()
            if quantification_init.weights:
                quant_lsv = lsvs[lidx]
                lsv_type = quant_lsv.type
                weights = np.tile(quantification_init.weights, (len(list_of_lsv), 1))

                psi = bootstrap_samples_calculation(quant_lsv, n_replica=num_exp,
                                                    name=quantification_init.names,
                                                    outdir=quantification_init.output,
                                                    nbins=quantification_init.nbins,
                                                    store_bootsamples=quantification_init.boots,
                                                    lock_array=quantification_init.lock_per_file,
                                                    fitfunc_r=fitfunc,
                                                    m_samples=quantification_init.m,
                                                    k_positions=quantification_init.k,
                                                    discardzeros=quantification_init.discardzeros,
                                                    trimborder=quantification_init.trimborder,
                                                    debug=quantification_init.debug)
                psi = np.array(psi) * weights
            else:
                lsv_id = list_of_lsv[lidx]
                psi, lsv_type = majiq_io.load_bootstrap_samples(list_of_lsv[lidx], f_list)
                psi = np.array(psi)

            if quantification_init.boots:
                continue

            num_ways = psi.shape[1]
            alpha_prior, beta_prior = get_prior_params(lsv_type, num_ways)
            post_psi = []
            mu_psi = []
            for p_idx in xrange(int(num_ways)):
                posterior = np.zeros(shape=nbins, dtype=np.float)
                mu_psi_m = []
                for m in xrange(quantification_init.m):

                    # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                    junc = np.array([psi[xx][p_idx][m] for xx in xrange(num_exp)])

                    all_sample = np.array([psi[xx][yy][m].sum() for xx in xrange(num_exp) for yy in xrange(num_ways)])
                    mu_psi_m.append(float(junc.sum()) / all_sample.sum())
                    data_given_psi = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                       alpha_prior[p_idx], beta_prior[p_idx]))
                    # normalizing
                    posterior += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))

                mu_psi.append(np.median(mu_psi_m))
                post_psi.append(posterior / quantification_init.m)
                if num_ways == 2:
                    break

            qm = QueueMessage(QUEUE_MESSAGE_PSI_RESULT, (post_psi, mu_psi, lsv_id), chnk)
            quantification_init.queue.put(qm, block=True)

        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        quantification_init.queue.put(qm, block=True)
        quantification_init.lock[chnk].acquire()
        quantification_init.lock[chnk].release()

        [xx.close() for xx in f_list]

    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()
##

conf = collections.namedtuple('conf', 'name output_dir silent debug markstacks')


class CalcPsi(BasicPipeline):

    def run(self):
        self.calcpsi()

    def calcpsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """
        majiq_utils.create_if_not_exists(self.logger_path)

        self.logger = majiq_utils.get_logger("%s/psi_majiq.log" % self.logger_path, silent=self.silent,
                                             debug=self.debug)

        self.logger.info("")
        self.logger.info("Command: %s" % self)
        self.logger.info("Running Psi ...")
        self.logger.info("GROUP: %s" % self.files)
        self.nbins = 40

        file_locks = None
        if self.only_boots:
            file_locks = [mp.Lock() for xx in self.files]
        lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()

        lsv_dict, lsv_types, lsv_summarized, meta = majiq_io.extract_lsv_summary(self.files)

        list_of_lsv = majiq_filter.merge_files_hdf5(lsv_dict=lsv_dict, lsv_summarized=lsv_summarized,
                                                    minnonzero=self.minpos, min_reads=self.minreads,
                                                    percent=self.min_exp, logger=self.logger)
        lchnksize = max(len(list_of_lsv)/self.nthreads, 1) + 1
        weights = self.calc_weights(self.weights, self.files, list_of_lsv, lock_arr, lchnksize,
                                    q, self.name)

        if self.only_boots:
            majiq_io.create_bootstrap_file(self.files, self.outDir, self.name, m=self.m)

        if len(list_of_lsv) > 0:
            [xx.acquire() for xx in lock_arr]
            pool = mp.Pool(processes=self.nthreads, initializer=quantification_init,
                           initargs=[q, lock_arr, self.outDir, self.name, self.silent, self.debug, self.nbins, self.m,
                                     self.k, self.discardzeros, self.trimborder, self.files, self.only_boots, weights,
                                     file_locks],
                           maxtasksperchild=1)

            pool.map_async(psi_quantification, majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            with Voila(get_quantifier_voila_filename(self.outDir, self.name), 'w') as out_h5p:
                out_h5p.add_metainfo(meta['genome'], self.name, meta['experiments'])
                in_h5p = h5py.File(self.files[0], 'r')
                queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=self.nthreads, logger=self.logger)
                in_h5p.close()

            pool.join()

        if self.only_boots:
            majiq_io.close_bootstrap_file(self.files, self.outDir, self.name, m=self.m)

        self.logger.info("PSI calculation for %s ended succesfully! "
                         "Result can be found at %s" % (self.name, self.outDir))
        self.logger.info("Alakazam! Done.")

