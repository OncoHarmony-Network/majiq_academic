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
from majiq.src.psi import prob_data_sample_given_psi, get_prior_params
from voila.io_voila import Voila

import collections


################################
# PSI calculation pipeline     #
################################


def boots_write(hg_grp, vals, lsv_idx):

    njunc = vals['samples'].shape[0]
    if lsv_idx + njunc > 2:
        shp = hg_grp['junctions'].shape
        shp_new = shp[0] + 5000
        hg_grp['junctions'].resize((shp_new, shp[1]))

    hg_grp['junctions'][lsv_idx:lsv_idx+njunc] = vals['samples']

    h_lsv = hg_grp.create_group("LSVs/%s" % vals['id'])
    h_lsv.attrs['id'] = vals['id']
    h_lsv.attrs['type'] = vals['type']
    h_lsv.attrs['coverage'] = hg_grp['junctions'].regionref[lsv_idx:lsv_idx + njunc]

    return lsv_idx+njunc


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
        num_exp = len(quantification_init.files)

        f_list = []
        # for eidx in np.arange(num_exp):
        #
        #     f_list.append(h5py.File(get_quantifier_norm_temp_files(quantification_init.output,
        #                                                           quantification_init.names, eidx)))

        for fname in quantification_init.files:
            f_list.append(h5py.File(fname))

        for lidx, lsv_id in enumerate(list_of_lsv):
            if lidx % 50 == 0:
                print "Event %d ..." % lidx
                sys.stdout.flush()

            lsv_samples = []
            for ff in f_list:
                lsvobj = ff['LSVs/%s' % lsv_id]
                lsv = ff[JUNCTIONS_DATASET_NAME][lsvobj.attrs['coverage']]
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
                    vals = {'samples': lsv_samples[ii], 'id': lsvobj.attrs['id'], 'type': lsvobj.attrs['type']}
                    file_name = '%s/%s.%d.boots.hdf5' % (quantification_init.output, quantification_init.names, ii)
                    quantification_init.lock_per_file[ii].acquire()
                    f = h5py.File(file_name, 'r+', compression='gzip', compression_opts=9)
                    lsv_idx = f.attrs['lsv_idx']
                    lsv_idx = boots_write(f, vals, lsv_idx)
                    f.attrs['lsv_idx'] = lsv_idx
                    f.close()
                    quantification_init.lock_per_file[ii].release()

                continue

            psi = np.array(lsv_samples)
            num_ways = psi.shape[1]
            alpha_prior, beta_prior = get_prior_params(lsvobj.attrs['type'], num_ways)
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


# def parse_and_norm_majiq(fname, replica_num, config):
#     try:
#         logger = majiq_utils.get_logger("%s/%s_%s.majiq.log" % (config.output_dir, config.name, replica_num),
#                                         silent=config.silent, debug=config.debug)
#         meta_info, lsv_junc = majiq_io.load_data_lsv(fname, config.name, logger)
#         filtered_lsv = majiq_norm.mark_stacks(lsv_junc, meta_info['fitfunc'], config.markstacks, logger)
#
#         sys.stdout.flush()
#
#         nlsvs = len(filtered_lsv[0])
#         effective_readlen = filtered_lsv[0][0].shape[1]
#         # effective_readlen = meta_info['effective_readlen']
#
#         with h5py.File(get_quantifier_norm_temp_files(config.output_dir, config.name, replica_num),
#                        'w', compression='gzip', compression_opts=9) as f:
#             f.create_dataset(JUNCTIONS_DATASET_NAME,
#                              (nlsvs*2, effective_readlen),
#                              maxshape=(None, effective_readlen))
#             f.attrs['meta_info'] = meta_info['group']
#             # f.attrs['effective_readlen'] = effective_readlen
#             f.attrs['fitfunc'] = meta_info['fitfunc']
#
#             old_shape = nlsvs
#             lsv_idx = 0
#
#             for lidx, lsv in enumerate(filtered_lsv[0]):
#                 h_lsv = f.create_group('LSV/%s' % filtered_lsv[1][lidx][0])
#                 h_lsv.attrs['id'] = filtered_lsv[1][lidx][0]
#                 h_lsv.attrs['type'] = filtered_lsv[1][lidx][1]
#                 filtered_lsv[1][lidx][2].copy(filtered_lsv[1][lidx][2], h_lsv, name='visual')
#
#                 all_vals = filtered_lsv[0][lidx]
#                 njunc = filtered_lsv[0][lidx].shape[0]
#
#                 if lsv_idx + njunc > old_shape:
#                     shp = f[JUNCTIONS_DATASET_NAME].shape
#                     shp_new = shp[0] + nlsvs
#                     old_shape = shp_new
#                     f[JUNCTIONS_DATASET_NAME].resize((shp_new, effective_readlen))
#
#                 f[JUNCTIONS_DATASET_NAME][lsv_idx:lsv_idx+njunc] = np.reshape(np.array(all_vals),
#                                                                                   (njunc, effective_readlen))
#                 h_lsv.attrs['coverage'] = f[JUNCTIONS_DATASET_NAME].regionref[lsv_idx:lsv_idx + njunc]
#                 lsv_idx += njunc
#
#             f[JUNCTIONS_DATASET_NAME].resize((lsv_idx, effective_readlen))
#     except Exception as e:
#         traceback.print_exc()
#         sys.stdout.flush()
#         raise


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

        # logger = self.logger
        # self.logger = None
        # pool = mp.Pool(processes=self.nthreads, maxtasksperchild=1)
        #
        # for fidx, fname in enumerate(self.files):
        #     if self.nthreads > 1:
        #         pool.apply_async(parse_and_norm_majiq, [fname, fidx, conf(output_dir=self.output, name=self.name,
        #                                                                   debug=self.debug, silent=self.silent,
        #                                                                   markstacks=self.markstacks)])
        #     else:
        #         parse_and_norm_majiq(fname, fidx, conf(output_dir=self.output, name=self.name,
        #                                                debug=self.debug, silent=self.silent,
        #                                                markstacks=self.markstacks))
        # pool.close()
        # pool.join()
        # self.logger = logger

        #list_of_lsv = majiq_filter.merge_files_hdf5([get_quantifier_norm_temp_files(self.output, self.name, xx)
        list_of_lsv = majiq_filter.merge_files_hdf5(self.files, self.minpos, self.minreads,
                                                    logger=self.logger)
        file_locks = None
        if self.only_boots:
            import datetime
            for ii, ff in enumerate(self.files):
                f = h5py.File('%s/%s.%d.boots.hdf5' % (self.output, self.name, ii),
                              'w', compression='gzip', compression_opts=9)
                f.create_dataset('junctions', (5000, self.m), maxshape=(None, self.m))
                # fill meta info
                f.attrs['sample_id'] = ff
                f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
                f.attrs['VERSION'] = VERSION
                f.attrs['lsv_idx'] = 0
                f.close()

            file_locks = [mp.Lock() for xx in self.files]
        lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()

        pool = mp.Pool(processes=self.nthreads, initializer=quantification_init,
                       initargs=[q, lock_arr, self.output, self.name, self.silent, self.debug, self.nbins, self.m,
                                 self.k, self.discardzeros, self.trimborder, self.files, self.only_boots,
                                 file_locks],
                       maxtasksperchild=1)

        lchnksize = max(len(list_of_lsv)/self.nthreads, 1) + 1
        [xx.acquire() for xx in lock_arr]

        if len(list_of_lsv) > 0:
            # psi_quantification((list_of_lsv, 0))
            pool.map_async(psi_quantification, majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            meta = majiq_io.read_meta_info(self.files)
            with Voila(get_quantifier_voila_filename(self.output, self.name), 'w') as out_h5p:
                out_h5p.add_metainfo(meta['genome'], self.name, meta['experiments'])
                in_h5p = h5py.File(self.files[0])
                queue_manager(in_h5p, out_h5p, lock_arr, q, num_chunks=self.nthreads, logger=self.logger)
                in_h5p.close()

            pool.join()

        if self.only_boots:
            for ii, ff in enumerate(self.files):
                file_name = '%s/%s.%d.boots.hdf5' % (self.output, self.name, ii)
                f = h5py.File(file_name, 'r+', compression='gzip', compression_opts=9)
                lsv_idx = f.attrs['lsv_idx']
                f['junctions'].resize((lsv_idx, self.m))
                f.close()

        self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s" % (self.name, self.output))
        self.logger.info("Alakazam! Done.")

