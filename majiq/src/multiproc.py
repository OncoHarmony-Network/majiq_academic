import Queue
import os
import sys
from multiprocessing import current_process

import majiq.src.io as majiq_io
from majiq.src import io as majiq_io
from majiq.src.constants import *
import majiq.src.utils as majiq_utils
import majiq.src.config as majiq_config

def parallel_lsv_child_calculation(func, args, tempdir, name, chunk, store=True):
    # try:
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)
    thread_logger = majiq_utils.get_logger("%s/majiq.%s.log" % (tempdir, chunk), silent=False)
    thread_logger.info("[Th %s]: START child,%s" % (chunk, current_process().name))

    args.append(thread_logger)
    results = func(*args)

    sys.stdout.flush()
    if store:
        thread_logger.info("[Th %s]: Saving ...%s " % (chunk, func.__name__))
        majiq_io.dump_bin_file(results, "%s/%s_th%s.%s.pickle" % (tempdir, name, chunk, func.__name__))


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
                        discardzeros, trimborder, num_exp, only_boots):

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
    quantification_init.num_exp = num_exp
    quantification_init.only_boots = only_boots


def queue_manager(input_h5dfp, output_h5dfp, lock_array, result_queue, num_chunks, meta_info=None, logger=None):

    nthr_count = 0
    posterior_matrix = []
    psi1 = []
    psi2 = []
    names = []

    lsv_idx = [0] * majiq_config.num_experiments
    while True:
        try:
            val = result_queue.get(block=True, timeout=10)
            if val.get_type() == QUEUE_MESSAGE_BUILD_LSV:
                for jdx, exp_idx in enumerate(majiq_config.tissue_repl[val.get_value()[1]]):
                    lsvobj = val.get_value()[0]
                    lsv_idx[exp_idx] = lsvobj.to_hdf5(hdf5grp=output_h5dfp[exp_idx],
                                                      lsv_idx=lsv_idx[exp_idx],
                                                      exp_idx=jdx)

            elif val.get_type() == QUEUE_MESSAGE_PSI_RESULT:
                posterior_matrix.append(val.get_value()[0])
                names.append([majiq_io.load_lsvgraphic_from_majiq(input_h5dfp, val.get_value()[-1])])

            elif val.get_type() == QUEUE_MESSAGE_DELTAPSI_RESULT:
                posterior_matrix.append(val.get_value()[0])
                psi1.append(val.get_value()[1])
                psi2.append(val.get_value()[2])
                names.append([majiq_io.load_lsvgraphic_from_majiq(input_h5dfp, val.get_value()[-1])])
                pass

            elif val.get_type() == QUEUE_MESSAGE_END_WORKER:
                lock_array[val.get_chunk()].release()
                nthr_count += 1
                if nthr_count >= num_chunks:
                    break

        except Queue.Empty:
            if nthr_count < num_chunks:
                continue
            break
