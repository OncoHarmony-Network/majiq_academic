import queue
import os
import sys

import multiprocessing as mp
from majiq.src import io as majiq_io
import majiq.src.logger as majiq_logger
from majiq.src.constants import *
from voila.vlsv import VoilaLsv
import traceback


def process_wrapper(args_vals):

    try:
        vals, chnk = args_vals
        logger = majiq_logger.get_logger("%s/%s.majiq.log" % (process_conf.outDir, chnk),
                                        silent=process_conf.silent, debug=process_conf.debug)

        process_conf.func(vals, chnk, process_conf, logger=logger)

    except:
        # majiq_utils.monitor('CHILD %s:: EXCEPT' % chnk)
        traceback.print_exc()
        sys.stdout.flush()
        raise

    finally:
        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        process_conf.queue.put(qm, block=True)
        process_conf.lock[chnk].acquire()
        process_conf.lock[chnk].release()
        process_conf.queue.close()
        majiq_logger.close_logger(logger)


def parallel_lsv_child_calculation(func, args, tempdir, name, chunk):
    # try:
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)
    thread_logger = majiq_logger.get_logger("%s/majiq.%s.log" % (tempdir, chunk), silent=False)
    thread_logger.info("[Th %s]: START child,%s" % (chunk, mp.current_process().name))

    args.append(thread_logger)
    func(*args)
    sys.stdout.flush()

def chunks(l, n_in, extra):
    """Yield successive n-sized chunks from l.
    :param l: list to be split
    :param n: max length of chunks
    """
    try:
        idx = -1
        n = int(n_in)
        for i in range(0, len(l), n):
            idx += 1
            if extra is not None:
                yield (l[i:i+n], extra[idx])
            else:
                yield l[i:i+n]
    except:
        print("ERROR: extra value has incorrect size %s" % idx, extra)
        raise


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


def process_conf(func, pipeline):
    process_conf.__dict__.update(pipeline.__dict__)
    process_conf.func = func


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
            sys.stdout.flush()
            if val.get_type() == QUEUE_MESSAGE_BUILD_JUNCTION:
                sys.stdout.flush()
                # majiq_io.add_lsv_to_bootstrapfile(output_h5dfp[val.get_value()[1]], val.get_value()[0])

            elif val.get_type() == QUEUE_MESSAGE_PSI_RESULT:

                lsv_graph = list_of_lsv_graphics[val.get_value()[-1]]
                output_h5dfp.add_lsv(VoilaLsv(bins_list=val.get_value()[0], means_psi1=val.get_value()[1],
                                              lsv_graphic=lsv_graph))

            elif val.get_type() == QUEUE_MESSAGE_DELTAPSI_RESULT:

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
            del val

        except queue.Empty:
            if nthr_count < num_chunks:
                continue
            break

