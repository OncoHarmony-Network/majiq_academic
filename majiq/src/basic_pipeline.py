import majiq.src.io as majiq_io
# from majiq.src.polyfitnb import fit_nb
import abc
import majiq.src.logger as majiq_logger
from majiq.src.psi import divs_from_bootsamples, calc_rho_from_divs, calc_local_weights
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, chunks, process_wrapper
from majiq.src.constants import *

import multiprocessing as mp


# ###############################
# Data loading and Boilerplate #
################################


# def bootstrap_samples_with_divs(list_of_lsv, chnk, conf, logger):
#
#     lsvs_to_work = majiq_io.get_extract_lsv_list(list_of_lsv, conf.file_list, aggr=False)
#     divs = divs_from_bootsamples(lsvs_to_work, n_replica=len(conf.file_list), pnorm=1, nbins=conf.nbins,
#                                  lsv_types=conf.lsv_type_dict)
#     qm = QueueMessage(QUEUE_MESSAGE_BOOTSTRAP, ([xx.id for xx in lsvs_to_work], divs), chnk)
#     conf.queue.put(qm, block=True)

def pipeline_run(pipeline):
    """ Exception catching for all the pipelines """
    try:
        return pipeline.run()
    except KeyboardInterrupt:
        if pipeline.logger:
            pipeline.logger.info("MAJIQ manually interrupted. Avada kedavra...")


class BasicPipeline:

    def __init__(self, args):
        """Basic configuration shared by all pipelines"""

        self.__dict__.update(args.__dict__)

        # if self.plotpath:
        #     majiq_logger.create_if_not_exists(self.plotpath)
        # self.logger_path = self.logger
        # if not self.logger_path:
        #     self.logger_path = self.outDir
        self.nthreads = args.nproc
        self.dpsi = False
        #self.psi_paths = []
        # try:
        #     self.replica_len = [len(self.files1), len(self.files2)]
        # except AttributeError:
        #     pass

    # def fitfunc(self, const_junctions):
    #     """Calculate the Negative Binomial function to sample from using the Constitutive events"""
    #     if self.debug:
    #         if self.logger is not None:
    #             self.logger.debug("Skipping fitfunc because --debug!")
    #         return np.poly1d([1, 0])
    #     else:
    #         if self.logger is not None:
    #             self.logger.debug("Fitting NB function with constitutive events...")
    #         return fit_nb(const_junctions, "%s/nbfit" % self.outDir, self.plotpath, logger=self.logger)

    # def calc_weights(self, weight_type, list_of_lsv, name, file_list, logger=None):
    #     self.file_list = file_list
    #     if weight_type.lower() == WEIGTHS_AUTO and len(self.file_list) >= 3:
    #         """ Calculate bootstraps samples and weights """
    #         nthreads = min(self.nthreads, len(list_of_lsv))
    #
    #         pool = mp.Pool(processes=nthreads, initializer=process_conf,
    #                        initargs=[bootstrap_samples_with_divs, self],
    #                        maxtasksperchild=1)
    #
    #         [xx.acquire() for xx in self.lock]
    #         pool.map_async(process_wrapper, chunks(list_of_lsv, nthreads))
    #         pool.close()
    #         divs = []
    #         lsvs = []
    #         queue_manager(output_h5dfp=None, lock_array=self.lock, result_queue=self.queue,
    #                       num_chunks=nthreads, out_inplace=(lsvs, divs), logger=logger)
    #         pool.join()
    #         # lsvs = {xx: vv for xx, vv in enumerate(lsvs)}
    #         divs = np.array(divs)
    #         rho = calc_rho_from_divs(divs, thresh=self.weights_threshold, alpha=self.weights_alpha,
    #                                  nreps=len(self.file_list), logger=self.logger)
    #
    #         wgts = calc_local_weights(divs, rho, self.local)
    #         if logger is not None:
    #             logger.info('Global weights for %s are: %s' % (name, ','.join([str(x) for x in rho])))
    #
    #         majiq_io.store_weights(lsvs, wgts, self.outDir, name)
    #         wgts = None
    #
    #     elif weight_type.lower() == WEIGTHS_NONE or len(self.file_list) < 3:
    #         wgts = np.ones(shape=(len(self.file_list)))
    #
    #     else:
    #         wgts = np.array([float(xx) for xx in weight_type.split(',')])
    #         if len(wgts) != len(self.file_list):
    #             logger.error('weights wrong arguments number of weights values is different than number of '
    #                          'specified replicas for that group')
    #     del self.file_list
    #     return wgts



    @abc.abstractmethod
    def run(self):
        """This is the entry point for all pipelines"""
        return