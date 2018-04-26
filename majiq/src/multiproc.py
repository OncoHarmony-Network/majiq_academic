import multiprocessing as mp
import os
import queue
import sys
import traceback

import psutil

import majiq.src.logger as majiq_logger
from majiq.src.constants import *

def process_wrapper(args_vals):
    try:
        vals, chnk = args_vals
        logger = majiq_logger.get_logger("%s/%s.majiq.log" % (process_conf.outDir, chnk),
                                         silent=process_conf.silent, debug=process_conf.debug)

        process_conf.func(vals, chnk, process_conf, logger=logger)
        logger.info('Finishing child, %s' % chnk)
        status = 0

    except Exception as e:
        logger.exception("Exception ocurred on %s" % process_conf.func.__name__)
        traceback.print_exc()
        sys.stdout.flush()
        status = -1
        raise

    finally:
        if process_conf.mem_profile:
            mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
            logger.info("Max Memory used %.2f MB" % mem_allocated)

        if process_conf.queue is not None:
            qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, status, chnk)
            logger.debug('SENDING END MESSAGE')
            process_conf.queue.put(qm, block=True)
            process_conf.lock[chnk].acquire()
            logger.debug('SENDING LOCK RELEASED')
            process_conf.lock[chnk].release()
            process_conf.queue.close()
            majiq_logger.close_logger(logger)


def parallel_lsv_child_calculation(func, args, tempdir, chunk):
    try:
        if not os.path.isdir(tempdir):
            os.mkdir(tempdir)
        thread_logger = majiq_logger.get_logger("%s/majiq.%s.log" % (tempdir, chunk), silent=False)
        thread_logger.info("[Th %s]: START child,%s" % (chunk, mp.current_process().name))

        args.append(thread_logger)
        func(*args)
        sys.stdout.flush()

    except Exception as e:
        thread_logger.exception("Exception ocurred on %s" % process_conf.func.__name__)
        traceback.print_exc()
        sys.stdout.flush()
        raise

    finally:
        majiq_logger.close_logger(thread_logger)


def chunks(l, n_chunks):
    """Yield successive n-sized chunks from l.
    :param l: list to be split
    :param n_chunks: total number of chunks
    """

    rem_len = len(l)
    n = 0
    prev_n = 0
    for ii in range(n_chunks):
        prev_n += n
        rem_len = rem_len - n
        n = int(rem_len / (n_chunks - ii))
        yield (l[prev_n:prev_n + n], ii)


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


def queue_manager(output, lock_array, result_queue, num_chunks, func, logger=None, **kwargs):
    nthr_count = 0
    kwargs['found'] = {}
    kwargs['gen_dict'] = {}
    err = False

    logger.debug('Start Queue Manager')
    while True:
        try:

            val = result_queue.get(block=True, timeout=10)
            if val.get_type() == QUEUE_MESSAGE_END_WORKER:
                status = val.get_value()
                logger.debug('WORKER DEATH MESSAGE RECEIVED %s (%s/%s)' % (val.get_chunk(), nthr_count, num_chunks))
                lock_array[val.get_chunk()].release()
                nthr_count += 1
                err = err or (status < 0)
                if nthr_count >= num_chunks:
                    break

            else:
                func(output=output, results=val.get_value(), msg_type=val.get_type(), extra=kwargs)

            del val

        except queue.Empty:
            if nthr_count < num_chunks:
                continue
            break

    return err
