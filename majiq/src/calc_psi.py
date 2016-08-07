import Queue
import multiprocessing as mp
import sys
import traceback
import h5py
import numpy as np
import scipy.misc

from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
import majiq.src.utils as majiq_utils
import majiq.src.filter as majiq_filter
import majiq.src.normalize as majiq_norm
import majiq.src.io as majiq_io
import majiq.src.sample as majiq_sample
from majiq.src.utils import chunks
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params
from majiq.src.constants import *


################################
# PSI calculation pipeline     #
################################


class QueueMessage:

    def __init__(self, msg_type, value, chunk):
        self.type = msg_type
        self.value = value
        self.chunk = chunk

    def is_closing(self):
        return self.type == -1

    def get_value(self):
        return self.value

    def get_chunk(self):
        return self.chunk

    def get_type(self):
        return self.type


def calcpsi(args):
    return pipeline_run(CalcPsi(args))


def psi_init(q, output, silent, debug, nbins, m, k, discardzeros, trimborder):

    psi_quantification.queue = q
    psi_quantification.output = output
    psi_quantification.silent = silent
    psi_quantification.debug = debug
    psi_quantification.nbins = nbins
    psi_quantification.m = m
    psi_quantification.k = k
    psi_quantification.discardzeros = discardzeros
    psi_quantification.trimborder = trimborder


def psi_quantification(args_vals):
    try:
        print "INIT MSG, ", args_vals
        created = mp.Process()
        chnk = created._identity[0]
        list_of_lsv, lock = args_vals
        logger = majiq_utils.get_logger("%s/%s.kk.majiq.log" % (psi_quantification.output, chnk),
                                        silent=psi_quantification.silent, debug=psi_quantification.debug)


        logger.info("Quantifying LSVs PSI.. %s" % chnk)
        nbins = psi_quantification.nbins

        fname = "%s/tmp.filtered.hdf5" % psi_quantification.output
        f = h5py.File(fname)
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
                                                                           m=psi_quantification.m,
                                                                           k=psi_quantification.k,
                                                                           discardzeros=psi_quantification.discardzeros,
                                                                           trimborder=psi_quantification.trimborder,
                                                                           fitted_one_over_r=f.attrs['fitfunc'][eidx],
                                                                           debug=psi_quantification.debug)

                lsv_samples.append(s_lsv)

            psi = np.array(lsv_samples)
            num_ways = len(lsv_samples[lidx])
            alpha_prior, beta_prior = get_prior_params(lsvobj.attrs['type'], num_ways)
            post_psi = []

            for p_idx in xrange(int(num_ways)):
                posterior = np.zeros(shape=nbins, dtype=np.float)

                for m in xrange(psi_quantification.m):
                    # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                    junc = np.array([psi[xx][p_idx][m] for xx in xrange(num_exp)])
                    all_sample = np.array([psi[xx][yy][m].sum() for xx in xrange(num_exp) for yy in xrange(num_ways)])
                    data_given_psi = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                       alpha_prior[p_idx], beta_prior[p_idx]))
                    # normalizing
                    posterior += np.exp(data_given_psi - scipy.misc.logsumexp(data_given_psi))

                post_psi.append(posterior / psi_quantification.m)
                if num_ways == 2:
                    break
            print "SEND LSV MSG ", chnk
            sys.stdout.flush()
            qm = QueueMessage(0, (post_psi, None), chnk)
            psi_quantification.queue.put(qm, block=True)

        print "SEND END MSG, ", lock
        sys.stdout.flush()
        qm = QueueMessage(-1, None, chnk)
        psi_quantification.queue.put(qm, block=True)
        lock[chnk].acquire()
        lock[chnk].release()

    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()




class CalcPsi(BasicPipeline):

    def run(self):
        self.calcpsi()

    def pre_psi(self, filename, logger=None):

        if logger is None:
            logger = majiq_utils.get_logger("%s/old_majiq.log" % self.output, silent=False)

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

    def queue_manage(self, lock_array, result_queue, first_id=0, meta_info=None, logger=None):

        fname = "%s/%s.psi.voila" % (self.output, self.name)

        # h5py.File(fname, 'w', compression='gzip', compression_opts=9)
        nthr_count = 0
        posterior_matrix = []
        names = []
        while True:
            try:
                val = result_queue.get(block=True, timeout=10)
                if val.get_type() == 0:
                    print "LSV RECEIVED"
                    sys.stdout.flush()
                    posterior_matrix.append(val.get_value()[0])
                    names.append(val.get_value()[1])

                elif val.get_type() == -1:
                    print "END WORKER"
                    sys.stdout.flush()

                    lock_array[val.get_chunk()].release()
                    nthr_count += 1

            except Queue.Empty:
                print "WAITING"
                sys.stdout.flush()
                if nthr_count < self.nthreads:
                    continue
                break

        pickle_path = "%s/%s_psigroup.pickle" % (self.output, self.name)
#        majiq_io.dump_lsvs_voila(pickle_path, posterior_matrix, names, meta_info)

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

        fname = "%s/tmp.filtered.hdf5" % self.output

        f = h5py.File(fname, 'w', compression='gzip', compression_opts=9)
        list_of_lsv = majiq_filter.quantifiable_in_group_to_hdf5(f, filtered_lsv, self.minpos, self.minreads,
                                                                 effective_readlen=61, logger=self.logger)

        f.attrs['num_exp'] = num_exp
        f.attrs['fitfunc'] = fitfunc
        f.close()

        lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()
        pool = mp.Pool(processes=self.nthreads, initializer=psi_init,
                       initargs=[q, self.output, self.silent, self.debug, self.nbins, self.m, self.k,
                                 self.discardzeros, self.trimborder], maxtasksperchild=1)

        lchnksize = max(len(list_of_lsv)/self.nthreads, 1)
        [xx.acquire() for xx in lock_arr]

        print "MAIN::", lock_arr
        if len(list_of_lsv) > 0:
            pool.map_async(psi_quantification, chunks(list_of_lsv, lchnksize, extra=lock_arr))
            pool.close()
            self.queue_manage(lock_arr, q, logger=self.logger)
            pool.join()

        self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s" % (self.name, self.output))
        self.logger.info("Alakazam! Done.")

