import h5py
import numpy as np
import scipy.misc
import sys
import multiprocessing as mp
import traceback

from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
import majiq.src.utils as majiq_utils
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
from majiq.src.io_utils import dump_bin_file, load_bin_file
import majiq.src.sample as majiq_sample
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params, gen_prior_matrix
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, quantification_init, queue_manager
from voila.io_voila import Voila
import collections


def deltapsi(args):
    return pipeline_run(DeltaPsi(args))


def deltapsi_quantification(args_vals):#, delta_prior_path, boots_sample=True, logger=None):
    try:

        list_of_lsv, chnk = args_vals
        logger = majiq_utils.get_logger("%s/%s.majiq.log" % (quantification_init.output, chnk),
                                        silent=quantification_init.silent, debug=quantification_init.debug)

        logger.info("Quantifying LSVs PSI.. %s" % chnk)
        num_exp = [len(quantification_init.files[0]), len(quantification_init.files[1])]
        f_list = [[], []]

        fitfunc = [[], []]
        for grp_idx in range(2):
            for fname in quantification_init.files[grp_idx]:

                f_list[grp_idx].append(majiq_io.open_hdf5_file(fname))
                fitfunc[grp_idx].append(f_list[grp_idx][-1].attrs['fitfunc'])
        prior_matrix = np.array(load_bin_file(get_prior_matrix_filename(quantification_init.output,
                                                                        quantification_init.names)))

        ones_n = np.ones(shape=(1, quantification_init.nbins), dtype=np.float)
        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print "Event %d ..." % lidx
                sys.stdout.flush()

            lsv_samples = [[], []]

            lsv_type = f_list[0][0]['LSVs/%s' % lsv_id].attrs['type']

            for grp_idx in range(2):
                for eidx, f in enumerate(f_list[grp_idx]):
                    lsvobj = f['LSVs/%s' % lsv_id]
                    assert lsv_type == lsvobj.attrs['type'], "LSV %s has different definition in groups" % lsv_id
                    lsv = f[JUNCTIONS_DATASET_NAME][lsvobj.attrs['coverage']]
                    m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(junction_list=lsv,
                                                                               m=quantification_init.m,
                                                                               k=quantification_init.k,
                                                                               discardzeros=quantification_init.discardzeros,
                                                                               trimborder=quantification_init.trimborder,
                                                                               fitted_one_over_r=fitfunc[grp_idx][eidx],
                                                                               debug=quantification_init.debug)
                    lsv_samples[grp_idx].append(s_lsv)

            psi1 = np.array(lsv_samples[0])
            psi2 = np.array(lsv_samples[1])

            prior_idx = 1 if 'i' in lsv_type else 0

            post_matrix = []
            posterior_psi1 = []
            posterior_psi2 = []
            num_ways = psi1.shape[1]
            mu_psi1 = []
            mu_psi2 = []
            alpha_prior, beta_prior = get_prior_params(lsv_type, num_ways)
            for p_idx in xrange(num_ways):

                posterior = np.zeros(shape=(quantification_init.nbins, quantification_init.nbins), dtype=np.float)
                post_psi1 = np.zeros(shape=quantification_init.nbins, dtype=np.float)
                post_psi2 = np.zeros(shape=quantification_init.nbins, dtype=np.float)
                mu_psi1_m = []
                mu_psi2_m = []
                for m in xrange(quantification_init.m):
                    # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                    junc = [psi1[xx][p_idx][m] for xx in xrange(num_exp[0])]
                    junc = np.array(junc)
                    all_sample = [psi1[xx][yy][m].sum() for xx in xrange(num_exp[0]) for yy in xrange(num_ways)]
                    all_sample = np.array(all_sample)
                    mu_psi1_m.append(float(junc.sum()) / all_sample.sum())
                    data_given_psi1 = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(),
                                                                        quantification_init.nbins,
                                                                        alpha_prior[p_idx], beta_prior[p_idx]))

                    psi_v1 = data_given_psi1.reshape(quantification_init.nbins, -1)
                    post_psi1 += np.exp(data_given_psi1 - scipy.misc.logsumexp(data_given_psi1))

                    junc = [psi2[xx][p_idx][m] for xx in xrange(num_exp[1])]
                    junc = np.array(junc)
                    all_sample = [psi2[xx][yy][m].sum() for xx in xrange(num_exp[1]) for yy in xrange(num_ways)]
                    all_sample = np.array(all_sample)
                    mu_psi2_m.append(float(junc.sum()) / all_sample.sum())
                    data_given_psi2 = np.log(prob_data_sample_given_psi(junc.sum(), all_sample.sum(),
                                                                        quantification_init.nbins,
                                                                        alpha_prior[p_idx], beta_prior[p_idx]))
                    post_psi2 += np.exp(data_given_psi2 - scipy.misc.logsumexp(data_given_psi2))
                    psi_v2 = data_given_psi2.reshape(-1, quantification_init.nbins)

                    A = (psi_v1 * ones_n + psi_v2 * ones_n.T) + np.log(prior_matrix[prior_idx])

                    posterior += np.exp(A - scipy.misc.logsumexp(A))
                mu_psi1.append(np.median(mu_psi1_m))
                mu_psi2.append(np.median(mu_psi2_m))
                post_matrix.append(posterior / quantification_init.m)
                posterior_psi1.append(post_psi1 / quantification_init.m)
                posterior_psi2.append(post_psi2 / quantification_init.m)
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

        # logger = self.logger
        # self.logger = None
        # pool = mp.Pool(processes=self.nthreads, maxtasksperchild=1)
        #
        # lsv_ids = []

        # for grp_idx, name in enumerate(self.names):
        #
        #     for fidx, fname in enumerate(files[grp_idx]):
        #         pool.apply_async(parse_and_norm_majiq, [fname, fidx, conf(output_dir=self.output, name=name,
        #                                                                   debug=self.debug, silent=self.silent,
        #                                                                   markstacks=self.markstacks)])
        # pool.close()
        # pool.join()
        # self.logger = logger

        list_of_lsv_group1 = majiq_filter.merge_files_hdf5(self.files1, self.minpos,
                                                           self.minreads, logger=self.logger)
        list_of_lsv_group2 = majiq_filter.merge_files_hdf5(self.files2, self.minpos,
                                                           self.minreads, logger=self.logger)

        list_of_lsv = list(set(list_of_lsv_group1).intersection(set(list_of_lsv_group2)))

        psi_space, prior_matrix = gen_prior_matrix(files, list_of_lsv_group1, list_of_lsv_group2, self.output,
                                                   prior_conf(self.iter, self.plotpath, self.breakiter, self.names,
                                                              self.binsize), numbins=self.nbins,
                                                   defaultprior=self.default_prior, logger=self.logger)

        self.logger.info("Saving prior matrix for %s..." % self.names)
        dump_bin_file(prior_matrix, get_prior_matrix_filename(self.output, self.names))

        lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()
        pool = mp.Pool(processes=self.nthreads, initializer=quantification_init,
                       initargs=[q, lock_arr, self.output, self.names, self.silent, self.debug, self.nbins, self.m,
                                 self.k, self.discardzeros, self.trimborder, files, False, None],
                       maxtasksperchild=1)
        lchnksize = max(len(list_of_lsv)/self.nthreads, 1) + 1
        [xx.acquire() for xx in lock_arr]

        if len(list_of_lsv) > 0:
            pool.map_async(deltapsi_quantification,
                           majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            meta1 = majiq_io.read_meta_info(self.files1)
            meta2 = majiq_io.read_meta_info(self.files2)
            with Voila(get_quantifier_voila_filename(self.output, self.names, deltapsi=True), 'w') as out_h5p:
                out_h5p.add_metainfo(meta1['genome'], group1=self.names[0], experiments1=meta1['experiments'],
                                     group2=self.names[1], experiments2=meta1['experiments'])

                in_h5p = h5py.File(files[0][0])
                queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=self.nthreads, logger=self.logger)
                in_h5p.close()


            pool.join()

        self.logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                          self.names[1],
                                                                                                          self.output))
        self.logger.info("Alakazam! Done.")
