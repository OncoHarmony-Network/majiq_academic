import multiprocessing as mp
import sys
import traceback

import h5py
import numpy as np
import scipy.misc

from majiq.src.polyfitnb import fit_nb
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.normalize as majiq_norm
import majiq.src.sample as majiq_sample
import majiq.src.utils as majiq_utils
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, quantification_init, queue_manager
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params
import collections


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
        if quantification_init.only_boots:
            logger.info(".. Only Bootstrap enabled .. %s" % chnk)

        nbins = quantification_init.nbins
        num_exp = quantification_init.num_exp

        f_list = []
        for eidx in np.arange(num_exp):
            f_list.append(h5py.File(get_quantifier_norm_temp_files(quantification_init.output,
                                                                   quantification_init.names, eidx)))
        if quantification_init.only_boots:
            res = [[] for xx in range(num_exp)]
            info_res = [[] for xx in range(num_exp)]

        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print "Event %d ..." % lidx
                sys.stdout.flush()

            lsv_samples = []
            for ff in f_list:
                lsvobj = ff['LSV/%s' % lsv_id]
                lsv = ff[LSV_JUNCTIONS_DATASET_NAME][lsvobj.attrs['coverage']]
                m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(junction_list=lsv,
                                                                           m=quantification_init.m,
                                                                           k=quantification_init.k,
                                                                           discardzeros=quantification_init.discardzeros,
                                                                           trimborder=quantification_init.trimborder,
                                                                           fitted_one_over_r=ff.attrs['fitfunc'],
                                                                           debug=quantification_init.debug)

                lsv_samples.append(s_lsv)

            if quantification_init.only_boots:
                for ii in range(num_exp):
                    res[ii].append(np.array(lsv_samples[ii]))
                    info_res[ii].append((lsvobj.attrs['id'], lsvobj.attrs['type']))
                continue

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

        if quantification_init.only_boots:
            import pickle
            for ii in range(num_exp):
                outfp = open('%s/%s.%d.boots.pickle' % (quantification_init.output, quantification_init.names, ii), 'w+')
                pickle.dump([res[ii], info_res[ii]], outfp)
                outfp.close()

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
def parse_and_norm_majiq(fname, replica_num, config):
    try:
        logger = majiq_utils.get_logger("%s/%s_%s.majiq.log" % (config.output_dir, config.name, replica_num),
                                        silent=config.silent, debug=config.debug)
        meta_info, lsv_junc = majiq_io.load_data_lsv(fname, config.name, logger)
        fitfunc = fit_nb(majiq_io.get_const_junctions(fname, logging=logger), "%s/nbfit" % config.output_dir,
                         None, logger=logger)
        filtered_lsv = majiq_norm.mark_stacks(lsv_junc, fitfunc, config.markstacks, logger)

        sys.stdout.flush()

        nlsvs = len(filtered_lsv[0])
        effective_readlen = filtered_lsv[0][0].shape[1]
        # effective_readlen = meta_info['effective_readlen']

        with h5py.File(get_quantifier_norm_temp_files(config.output_dir, config.name, replica_num),
                       'w', compression='gzip', compression_opts=9) as f:
            f.create_dataset(LSV_JUNCTIONS_DATASET_NAME,
                             (nlsvs*2, effective_readlen),
                             maxshape=(None, effective_readlen))
            f.attrs['meta_info'] = meta_info['group']
            # f.attrs['effective_readlen'] = effective_readlen
            f.attrs['fitfunc'] = fitfunc

            old_shape = nlsvs
            lsv_idx = 0

            for lidx, lsv in enumerate(filtered_lsv[0]):
                h_lsv = f.create_group('LSV/%s' % filtered_lsv[1][lidx][0])
                h_lsv.attrs['id'] = filtered_lsv[1][lidx][0]
                h_lsv.attrs['type'] = filtered_lsv[1][lidx][1]
                filtered_lsv[1][lidx][2].copy(filtered_lsv[1][lidx][2], h_lsv, name='visuals')

                all_vals = filtered_lsv[0][lidx]
                njunc = filtered_lsv[0][lidx].shape[0]

                if lsv_idx + njunc > old_shape:
                    shp = f[LSV_JUNCTIONS_DATASET_NAME].shape
                    shp_new = shp[0] + nlsvs
                    old_shape = shp_new
                    f[LSV_JUNCTIONS_DATASET_NAME].resize((shp_new, effective_readlen))

                f[LSV_JUNCTIONS_DATASET_NAME][lsv_idx:lsv_idx+njunc] = np.reshape(np.array(all_vals),
                                                                                  (njunc, effective_readlen))
                h_lsv.attrs['coverage'] = f[LSV_JUNCTIONS_DATASET_NAME].regionref[lsv_idx:lsv_idx + njunc]
                lsv_idx += njunc

            f[LSV_JUNCTIONS_DATASET_NAME].resize((lsv_idx, effective_readlen))
    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise

class CalcPsi(BasicPipeline):

    def run(self):
        self.calcpsi()

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

        logger = self.logger
        self.logger = None
        pool = mp.Pool(processes=self.nthreads, maxtasksperchild=1)

        for fidx, fname in enumerate(self.files):
            pool.apply_async(parse_and_norm_majiq, [fname, fidx, conf(output_dir=self.output, name=self.name,
                                                                           debug=self.debug, silent=self.silent,
                                                                           markstacks=self.markstacks)])

            # CalcPsi.parse_and_norm_majiq(fname, fidx, conf(output_dir=self.output, name=self.name,
            #                                                debug=self.debug, silent=self.silent,
            #                                                markstacks=self.markstacks,
            #                                                fitfunc=self.fitfunc))
        pool.close()
        pool.join()
        self.logger = logger

        list_of_lsv = majiq_filter.merge_files_hdf5([get_quantifier_norm_temp_files(self.output, self.name, xx)
                                                     for xx in xrange(len(self.files))], self.minpos, self.minreads,
                                                    logger=self.logger)

        lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()
        processes = self.nthreads if not self.only_boots else 1

        pool = mp.Pool(processes=processes, initializer=quantification_init,
                       initargs=[q, lock_arr, self.output, self.name, self.silent, self.debug, self.nbins, self.m,
                                 self.k, self.discardzeros, self.trimborder, len(self.files), self.only_boots],
                       maxtasksperchild=1)

        lchnksize = max(len(list_of_lsv)/processes, 1) + 1
        [xx.acquire() for xx in lock_arr]

        if len(list_of_lsv) > 0:
            pool.map_async(psi_quantification, majiq_utils.chunks(list_of_lsv, lchnksize, extra=range(processes)))
            pool.close()

            out_h5p = h5py.File(get_quantifier_voila_filename(self.output, self.name),
                                'w', compression='gzip', compression_opts=9)
            in_h5p = h5py.File(get_quantifier_norm_temp_files(self.output, self.name, 0))

            queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=processes, logger=self.logger)

            in_h5p.close()
            out_h5p.close()

            pool.join()

        self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s" % (self.name, self.output))
        self.logger.info("Alakazam! Done.")

