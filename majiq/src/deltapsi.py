import h5py
import numpy as np
import scipy.misc
import sys
import multiprocessing as mp
import Queue
import traceback

from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
import majiq.src.utils as majiq_utils
import majiq.src.filter as majiq_filter
import majiq.src.normalize as majiq_norm
import majiq.src.io as majiq_io
import majiq.src.sample as majiq_sample
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params, gen_prior_matrix
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, quantification_init, queue_manager

def deltapsi(args):
    return pipeline_run(DeltaPsi(args))


def deltapsi_quantification(args_vals):#, delta_prior_path, boots_sample=True, logger=None):
    try:

        list_of_lsv, chnk = args_vals
        logger = majiq_utils.get_logger("%s/%s.majiq.log" % (quantification_init.output, chnk),
                                        silent=quantification_init.silent, debug=quantification_init.debug)

        logger.info("Quantifying LSVs PSI.. %s" % chnk)
        filesp = [h5py.File(get_quantifier_temp_filename(quantification_init.output, quantification_init.names[0])),
                  h5py.File(get_quantifier_temp_filename(quantification_init.output, quantification_init.names[1]))]

        num_exp = [filesp[0].attrs['num_exp'], filesp[1].attrs['num_exp']]
        prior_matrix = np.array(majiq_io.load_bin_file(get_prior_matrix_filename(quantification_init.output,
                                                       quantification_init.names)))

        ones_n = np.ones(shape=(1, quantification_init.nbins), dtype=np.float)

        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print "Event %d ..." % lidx
                sys.stdout.flush()

            lsv_samples1 = []
            lsv_samples2 = []

            lsv_type = filesp[0][lsv_id].attrs['type']

            logger.info("Bootstrapping for all samples...")
            for grp_idx in range(2):
                f = filesp[grp_idx]
                lsvobj = f[lsv_id]

                assert lsv_type == lsvobj.attrs['type'], "LSV %s has different definition in groups" % lsv_id

                for eidx in np.arange(num_exp[grp_idx]):
                    lsv = f[grp_idx][LSV_JUNCTIONS_DATASET_NAME][lsvobj.attrs['coverage']][grp_idx, :, eidx, :]
                    m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(junction_list=lsv,
                                                                               m=quantification_init.m,
                                                                               k=quantification_init.k,
                                                                               discardzeros=quantification_init.discardzeros,
                                                                               trimborder=quantification_init.trimborder,
                                                                               fitted_one_over_r=f.attrs['fitfunc'][eidx],
                                                                               debug=quantification_init.debug)
                    if grp_idx == 0:
                        lsv_samples1.append(s_lsv)
                    else:
                        lsv_samples2.append(s_lsv)

            psi1 = np.array(lsv_samples1)
            psi2 = np.array(lsv_samples2)

            alpha_prior, beta_prior = get_prior_params(lsv_type, num_ways)

            if 'i' in lsv_type:
                prior_idx = 1
            else:
                prior_idx = 0

            post_matrix = []
            posterior_psi1 = []
            posterior_psi2 = []

            num_ways = len(psi1)

            for p_idx in xrange(num_ways):

                posterior = np.zeros(shape=(quantification_init.nbins, quantification_init.nbins), dtype=np.float)
                post_psi1 = np.zeros(shape=quantification_init.nbins, dtype=np.float)
                post_psi2 = np.zeros(shape=quantification_init.nbins, dtype=np.float)
                for m in xrange(quantification_init.m):
                    # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                    junc = [psi1[xx][p_idx][m] for xx in xrange(num_exp[0])]
                    junc = np.array(junc)
                    all_sample = [psi1[xx][yy][m].sum() for xx in xrange(num_exp[0]) for yy in xrange(num_ways)]
                    all_sample = np.array(all_sample)
                    data_given_psi1 = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(),
                                                                        quantification_init.nbins,
                                                                        alpha_prior[p_idx], beta_prior[p_idx]))

                    psi_v1 = data_given_psi1.reshape(quantification_init.nbins, -1)
                    post_psi1 += (data_given_psi1 - scipy.misc.logsumexp(data_given_psi1))

                    junc = [psi2[xx][p_idx][m] for xx in xrange(num_exp[1])]
                    junc = np.array(junc)
                    all_sample = [psi2[xx][yy][m].sum() for xx in xrange(num_exp[1]) for yy in xrange(num_ways)]
                    all_sample = np.array(all_sample)
                    data_given_psi2 = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(),
                                                                        quantification_init.nbins,
                                                                        alpha_prior[p_idx], beta_prior[p_idx]))
                    post_psi2 += (data_given_psi2 - scipy.misc.logsumexp(data_given_psi2))
                    psi_v2 = data_given_psi2.reshape(-1, quantification_init.nbins)

                    A = (psi_v1 * ones_n + psi_v2 * ones_n.T) + np.log(prior_matrix[prior_idx])

                    posterior += np.exp(A - scipy.misc.logsumexp(A))

                post_matrix.append(posterior / quantification_init.m)
                posterior_psi1.append(post_psi1 / quantification_init.m)
                posterior_psi2.append(post_psi2 / quantification_init.m)
                if num_ways == 2:
                    break

            qm = QueueMessage(QUEUE_MESSAGE_DELTAPSI_RESULT, (post_matrix, posterior_psi1, posterior_psi2, lsv_id), chnk)
            quantification_init.queue.put(qm, block=True)

        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        quantification_init.queue.put(qm, block=True)
        quantification_init.lock[chnk].acquire()
        quantification_init.lock[chnk].release()

    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


class DeltaPsi(BasicPipeline):

    def run(self):
        self.deltapsi()

    def deltapsi(self):
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
        self.nbins = 20

        lsv_ids = []
        for grp_idx, name in self.names:

            num_exp = len(self.files[grp_idx])
            meta_info = [0] * num_exp
            filtered_lsv = [None] * num_exp
            fitfunc = [None] * num_exp
            for ii, fname in enumerate(self.files[grp_idx]):
                meta_info[ii], lsv_junc = majiq_io.load_data_lsv(fname, name, self.logger)

                fitfunc[ii] = self.fitfunc(majiq_io.get_const_junctions(fname, logging=self.logger))
                filtered_lsv[ii] = majiq_norm.mark_stacks(lsv_junc, fitfunc[ii], self.markstacks, self.logger)

            f = h5py.File(get_quantifier_temp_filename(self.output, name), 'w', compression='gzip', compression_opts=9)
            lsv_ids.append(majiq_filter.quantifiable_in_group_to_hdf5(f, filtered_lsv, self.minpos, self.minreads,
                                                                      effective_readlen=61, logger=self.logger))

            f.attrs['num_exp'] = num_exp
            f.attrs['fitfunc'] = fitfunc
            f.close()

        list_of_lsv = list(lsv_ids[0].intersection(lsv_ids[1]))


        #group1, group2 = combine_for_priormatrix(matched_lsv[0], matched_lsv[1], matched_info, num_exp)
        psi_space, prior_matrix = gen_prior_matrix(self, group1, group2, self.output, numbins=self.nbins,
                                                   defaultprior=self.default_prior)

        self.logger.info("Saving prior matrix for %s..." % self.names)
        majiq_io.dump_bin_file(prior_matrix, get_prior_matrix_filename(self.output,
                                                                       self.names))



        lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()
        pool = mp.Pool(processes=self.nthreads, initializer=quantification_init,
                       initargs=[q, lock_arr, self.output, self.silent, self.debug, self.nbins, self.m, self.k,
                                 self.discardzeros, self.trimborder], maxtasksperchild=1)
        lchnksize = max(len(list_of_lsv)/self.nthreads, 1)
        [xx.acquire() for xx in lock_arr]

        if len(list_of_lsv) > 0:
            pool.map_async(deltapsi_quantification,
                           majiq_utils.chunks(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()

            out_h5p = h5py.File(get_quantifier_voila_filename(self.output, self.names, deltapsi=True),
                                'w', compression='gzip', compression_opts=9)
            in_h5p = h5py.File(get_quantifier_temp_filename(self.output, name))
            queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=self.nthreads, logger=self.logger)
            in_h5p.close()
            out_h5p.close()

            pool.join()

        self.logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                          self.names[1],
                                                                                                          self.output))
        self.logger.info("Alakazam! Done.")