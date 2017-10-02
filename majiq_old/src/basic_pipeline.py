import majiq.src.io as majiq_io
from majiq.src.polyfitnb import fit_nb
import abc

# from numpy.ma import masked_less
import majiq.src.utils as majiq_utils
from majiq.src.psi import divs_from_bootsamples, calc_rho_from_divs, calc_local_weights
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager
from majiq.src.constants import *
import sys
import traceback
import multiprocessing as mp


# ###############################
# Data loading and Boilerplate #
################################


# def get_clean_raw_reads(matched_info, matched_lsv, outdir, names, num_exp):
#     res = []
#     for eidx in range(num_exp):
#         for ldx, lsv in enumerate(matched_info):
#             jlist = masked_less(matched_lsv[ldx][eidx], 0)
#
#             num = jlist.sum(axis=1)
#             res.append([lsv[1], num.data])
#     majiq_io.dump_bin_file(res, '%s/clean_reads.%s.pkl' % (outdir, names))


def bootstrap_samples_with_divs(args_vals):

    try:
        list_of_lsv, chnk = args_vals
        logger = majiq_utils.get_logger("%s/%s.majiq.log" % (process_conf.output, chnk),
                                        silent=process_conf.silent, debug=process_conf.debug)

        lsvs_to_work, fitfunc_r = majiq_io.get_extract_lsv_list(list_of_lsv, process_conf.files)

        divs = divs_from_bootsamples(lsvs_to_work, fitfunc_r=fitfunc_r,
                                     n_replica=len(process_conf.files), pnorm=1, m_samples=process_conf.m,
                                     k_positions=process_conf.k, discardzeros=process_conf.discardzeros,
                                     trimborder=process_conf.trimborder, debug=False,
                                     nbins=process_conf.nbins, store_bootsamples=True,
                                     lock_array=process_conf.lock_per_file,
                                     outdir=process_conf.output, name=process_conf.names)

        qm = QueueMessage(QUEUE_MESSAGE_BOOTSTRAP, (lsvs_to_work, divs), chnk)
        process_conf.queue.put(qm, block=True)

        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        process_conf.queue.put(qm, block=True)
        process_conf.lock[chnk].acquire()
        process_conf.lock[chnk].release()

    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


def pipeline_run(pipeline):
    """ Exception catching for all the pipelines """
    try:
        return pipeline.run()
    except KeyboardInterrupt:
        if pipeline.logger:
            pipeline.logger.info("MAJIQ manually interrupted. Avada kedavra...")


class BasicPipeline:
    #
    # boots_conf = collections.namedtuple('boots_conf', 'm, k, discardzeros, trimborder')
    # basic_conf = collections.namedtuple('basic_conf', 'output, nbins, silent, debug')

    def __init__(self, args):
        """Basic configuration shared by all pipelines"""

        self.__dict__.update(args.__dict__)

        if self.plotpath:
            majiq_utils.create_if_not_exists(self.plotpath)
        self.logger_path = self.logger
        if not self.logger_path:
            self.logger_path = self.outDir

        self.nthreads = args.nthreads
        self.psi_paths = []
        self.nchunks = self.nthreads
        try:
            self.replica_len = [len(self.files1), len(self.files2)]
        except AttributeError:
            pass

    def fitfunc(self, const_junctions):
        """Calculate the Negative Binomial function to sample from using the Constitutive events"""
        if self.debug:
            if self.logger is not None:
                self.logger.debug("Skipping fitfunc because --debug!")
            return np.poly1d([1, 0])
        else:
            if self.logger is not None:
                self.logger.debug("Fitting NB function with constitutive events...")
            return fit_nb(const_junctions, "%s/nbfit" % self.outDir, self.plotpath, logger=self.logger)

    def calc_weights(self, weight_type, file_list, list_of_lsv, lock_arr, lchnksize, q, name, store=True, logger=None):

        if weight_type.lower() == WEIGTHS_AUTO and len(file_list) >= 3:
            """ Calculate bootstraps samples and weights """

            pool = mp.Pool(processes=self.nthreads, initializer=process_conf, initargs=[self, q, lock_arr, None],
                           maxtasksperchild=1)

            [xx.acquire() for xx in lock_arr]
            pool.map_async(bootstrap_samples_with_divs,
                           majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            divs = []
            lsvs = []
            queue_manager(input_h5dfp=None, output_h5dfp=None, lock_array=lock_arr, result_queue=q,
                          num_chunks=self.nthreads, out_inplace=(lsvs, divs),
                          logger=self.logger)
            pool.join()
            lsvs = {xx: vv for xx, vv in enumerate(lsvs)}
            divs = np.array(divs)
            rho = calc_rho_from_divs(divs, thresh=self.weights_threshold, alpha=self.weights_alpha,
                                     nreps=len(file_list), logger=self.logger)

            wgts = calc_local_weights(divs, rho, self.local)
            if store:

                majiq_io.store_weights_bootstrap(lsvs, wgts, file_list, self.outDir, name)
                wgts = None
            else:
                wgts = rho

        elif weight_type.lower() == WEIGTHS_NONE or len(file_list) < 3:
            wgts = np.ones(shape=(len(file_list)))

        else:
            wgts = np.array([float(xx) for xx in weight_type.split(',')])
            if len(wgts) != len(file_list):
                logger.error('weights wrong arguments number of weights values is different than number of '
                                  'specified replicas for that group')

        return wgts



    @abc.abstractmethod
    def run(self):
        """This is the entry point for all pipelines"""
        return