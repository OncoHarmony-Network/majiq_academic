import multiprocessing as mp
import sys
import traceback

import h5py
import numpy as np
import scipy.misc

import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.normalize as majiq_norm
import majiq.src.sample as majiq_sample
import majiq.src.utils as majiq_utils
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, quantification_init, queue_manager
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params


################################
# PSI calculation pipeline     #
################################


def calcpsi(args):
    return pipeline_run(CalcPsi(args))


def psi_quantification(args_vals):
    try:

        list_of_lsv, chnk = args_vals
        logger = majiq_utils.get_logger("%s/%s.majiq.log" % (quantification_init.output, chnk),
                                        silent=quantification_init.silent, debug=quantification_init.debug)

        logger.info("Quantifying LSVs PSI.. %s" % chnk)
        nbins = quantification_init.nbins

        f = h5py.File(get_quantifier_temp_filename(quantification_init.output, quantification_init.names))
        num_exp = f.attrs['num_exp']

        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print "Event %d ..." % lidx
                sys.stdout.flush()

            lsv_samples = []
            lsvobj = f[lsv_id]
            for eidx in np.arange(num_exp):
                lsv = f[LSV_JUNCTIONS_DATASET_NAME][lsvobj.attrs['coverage']][:, eidx, :]

                m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(junction_list=lsv,
                                                                           m=quantification_init.m,
                                                                           k=quantification_init.k,
                                                                           discardzeros=quantification_init.discardzeros,
                                                                           trimborder=quantification_init.trimborder,
                                                                           fitted_one_over_r=f.attrs['fitfunc'][eidx],
                                                                           debug=quantification_init.debug)

                lsv_samples.append(s_lsv)

            psi = np.array(lsv_samples)
            num_ways = len(lsv_samples)
            alpha_prior, beta_prior = get_prior_params(lsvobj.attrs['type'], num_ways)
            post_psi = []

            for p_idx in xrange(int(num_ways)):
                posterior = np.zeros(shape=nbins, dtype=np.float)

                for m in xrange(quantification_init.m):
                    # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                    junc = np.array([psi[xx][p_idx][m] for xx in xrange(num_exp)])
                    all_sample = np.array([psi[xx][yy][m].sum() for xx in xrange(num_exp) for yy in xrange(num_ways)])
                    data_given_psi = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                       alpha_prior[p_idx], beta_prior[p_idx]))
                    # normalizing
                    posterior += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))

                post_psi.append(posterior / quantification_init.m)
                if num_ways == 2:
                    break

            qm = QueueMessage(QUEUE_MESSAGE_PSI_RESULT, (post_psi, lsv_id), chnk)
            quantification_init.queue.put(qm, block=True)

        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        quantification_init.queue.put(qm, block=True)
        quantification_init.lock[chnk].acquire()
        quantification_init.lock[chnk].release()

    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


class CalcPsi(BasicPipeline):

    def run(self):
        self.calcpsi()

    def pre_psi(self, filename, logger=None):

        if logger is None:
            logger = majiq_utils.get_logger("%s/majiq.log" % self.output, silent=False)

        self.logger = logger
        num_exp = len(self.files)
        meta_info = [0] * num_exp
        filtered_lsv = [None] * num_exp
        fitfunc = [None] * num_exp
        for ii, fname in enumerate(self.files):
            meta_info[ii], lsv_junc = majiq_io.load_data_lsv(fname, self.name, logger)
            meta_info[ii], lsv_junc, const = majiq_io.get_const_junctions(fname, self.name, logger)
            fitfunc[ii] = self.fitfunc(const[0])
            filtered_lsv[ii] = majiq_norm.mark_stacks(lsv_junc, fitfunc[ii], self.markstacks, self.logger)

        matched_lsv, matched_info = majiq_filter.quantifiable_in_group(filtered_lsv, self.minpos, self.minreads,
                                                                       logger=logger)

    def calcpsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """

        self.logger = majiq_utils.get_logger("%s/psi_majiq.log" % self.logger_path, silent=self.silent,
                                             debug=self.debug)

        self.logger.info("")
        self.logger.info("Command: %s" % self)
        self.logger.info("Running Psi ...")
        self.logger.info("GROUP: %s" % self.files)
        self.nbins = 40

        num_exp = len(self.files)
        meta_info = [0] * num_exp
        filtered_lsv = [None] * num_exp
        fitfunc = [None] * num_exp
        for ii, fname in enumerate(self.files):
            meta_info[ii], lsv_junc = majiq_io.load_data_lsv(fname, self.name, self.logger)

            fitfunc[ii] = self.fitfunc(majiq_io.get_const_junctions(fname, logging=self.logger))
            filtered_lsv[ii] = majiq_norm.mark_stacks(lsv_junc, fitfunc[ii], self.markstacks, self.logger)

        f = h5py.File(get_quantifier_temp_filename(self.output, self.name),
                      'w', compression='gzip', compression_opts=9)
        list_of_lsv = majiq_filter.quantifiable_in_group_to_hdf5(f, filtered_lsv, self.minpos, self.minreads,
                                                                 effective_readlen=61, logger=self.logger)

        f.attrs['num_exp'] = num_exp
        f.attrs['fitfunc'] = fitfunc
        f.close()

        lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()
        pool = mp.Pool(processes=self.nthreads, initializer=quantification_init,
                       initargs=[q, lock_arr, self.output, self.name, self.silent, self.debug, self.nbins, self.m, self.k,
                                 self.discardzeros, self.trimborder], maxtasksperchild=1)
        lchnksize = max(len(list_of_lsv)/self.nthreads, 1)
        [xx.acquire() for xx in lock_arr]

        if len(list_of_lsv) > 0:
            pool.map_async(psi_quantification, chunks(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()

            out_h5p = h5py.File(get_quantifier_voila_filename(self.output, self.name),
                                'w', compression='gzip', compression_opts=9)
            in_h5p = h5py.File(get_quantifier_temp_filename(self.output, self.name))

            queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=self.nthreads, logger=self.logger)

            in_h5p.close()
            out_h5p.close()

            pool.join()

        self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s" % (self.name, self.output))
        self.logger.info("Alakazam! Done.")

