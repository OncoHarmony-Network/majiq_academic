import Queue
import os
import sys

import multiprocessing as mp
from majiq.src import io as majiq_io
from majiq.src.constants import *
import majiq.src.utils as majiq_utils
from majiq.src.config import Config
from voila.splice_graphics import LsvGraphic
from voila.vlsv import VoilaLsv


def parallel_lsv_child_calculation(func, args, tempdir, name, chunk, store=True):
    # try:
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)
    thread_logger = majiq_utils.get_logger("%s/majiq.%s.log" % (tempdir, chunk), silent=False)
    thread_logger.info("[Th %s]: START child,%s" % (chunk, mp.current_process().name))

    args.append(thread_logger)
    results = func(*args)

    sys.stdout.flush()
    if store:
        thread_logger.info("[Th %s]: Saving ...%s " % (chunk, func.__name__))
        majiq_io.dump_bin_file(results, "%s/%s_th%s.%s.pickle" % (tempdir, name, chunk, func.__name__))


# def pool_process(func, iter_args, nthreads, initializer, init_args,
#                  input_h5dfp=None, output_h5dfp=None, out_inplace=None, logger=None):
#     pool = mp.Pool(processes=nthreads, initializer=initializer,
#                    initargs=init_args,
#                    maxtasksperchild=1)
#     lchnksize = max(len(iter_args) / nthreads, 1) + 1
#     [xx.acquire() for xx in lock_arr]
#     pool.map_async(func, majiq_utils.chunks2(iter_args, lchnksize, extra=range(nthreads)))
#     pool.close()
#     queue_manager(input_h5dfp=input_h5dfp, output_h5dfp=output_h5dfp, lock_array=lock_arr, result_queue=q,
#                   num_chunks=nthreads, out_inplace=out_inplace,
#                   logger=logger)


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


def quantification_init(q, lock, output, names, silent, debug, nbins, m, k,
                        discardzeros, trimborder, files, only_boots, weights, lock_per_file):
    quantification_init.lock = lock
    quantification_init.queue = q
    quantification_init.output = output
    quantification_init.names = names
    quantification_init.silent = silent
    quantification_init.debug = debug
    quantification_init.nbins = nbins
    quantification_init.m = m
    quantification_init.k = k
    quantification_init.discardzeros = discardzeros
    quantification_init.trimborder = trimborder
    quantification_init.files = files
    quantification_init.boots = only_boots
    quantification_init.weights = weights
    quantification_init.lock_per_file = lock_per_file


def queue_manager(input_h5dfp, output_h5dfp, lock_array, result_queue, num_chunks, meta_info=None, num_exp=0,
                  out_inplace=None, logger=None, list_of_lsv_graphics={}):

    nthr_count = 0
    # posterior_matrix = []
    # psi1 = []
    # psi2 = []
    # names = []
    lsv_idx = [0] * num_exp
    while True:
        try:
            val = result_queue.get(block=True, timeout=10)
            if val.get_type() == QUEUE_MESSAGE_BUILD_LSV:
                majiq_config = Config()
                for jdx, exp_idx in enumerate(majiq_config.tissue_repl[val.get_value()[1]]):
                    lsvobj = val.get_value()[0]
                    lsv_idx[exp_idx] = lsvobj.to_hdf5(hdf5grp=output_h5dfp[exp_idx],
                                                      lsv_idx=lsv_idx[exp_idx],
                                                      exp_idx=jdx)

            elif val.get_type() == QUEUE_MESSAGE_PSI_RESULT:
                lsv_graph = list_of_lsv_graphics[val.get_value()[-1]]
                output_h5dfp.add_lsv(VoilaLsv(bins_list=val.get_value()[0], means_psi1=val.get_value()[1],
                                              lsv_graphic=lsv_graph))

            elif val.get_type() == QUEUE_MESSAGE_DELTAPSI_RESULT:

                # lsv_graph = LsvGraphic.easy_from_hdf5(majiq_io.load_lsvgraphic_from_majiq(input_h5dfp,
                #                                                                           val.get_value()[-1]))
                lsv_graph = list_of_lsv_graphics[val.get_value()[-1]]
                output_h5dfp.add_lsv(VoilaLsv(bins_list=val.get_value()[0], lsv_graphic=lsv_graph,
                                              psi1=val.get_value()[1], psi2=val.get_value()[2],
                                              means_psi1=val.get_value()[3], means_psi2=val.get_value()[4]))

            elif val.get_type() == QUEUE_MESSAGE_HETER_DELTAPSI:
                lsv_graph = list_of_lsv_graphics[val.get_value()[-1]]

                output_h5dfp.add_lsv(VoilaLsv(bins_list=None, lsv_graphic=lsv_graph, psi1=None, psi2=None,
                                              means_psi1=None, means_psi2=None, het=val.get_value()[0]))


            elif val.get_type() == QUEUE_MESSAGE_BOOTSTRAP:
                out_inplace[0].extend(val.get_value()[0])
                out_inplace[1].extend(val.get_value()[1])

            elif val.get_type() == QUEUE_MESSAGE_END_WORKER:
                lock_array[val.get_chunk()].release()
                nthr_count += 1
                if nthr_count >= num_chunks:
                    break

        except Queue.Empty:
            if nthr_count < num_chunks:
                continue
            break

