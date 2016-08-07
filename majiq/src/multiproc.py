import os
import sys
from multiprocessing import current_process

import old_majiq.src.io as majiq_io
from majiq.src.utils import get_logger


def parallel_lsv_child_calculation(func, args, tempdir, name, chunk, store=True):
    # try:
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)
    thread_logger = get_logger("%s/old_majiq.%s.log" % (tempdir, chunk), silent=False)
    thread_logger.info("[Th %s]: START child,%s" % (chunk, current_process().name))

    args.append(thread_logger)
    results = func(*args)

    sys.stdout.flush()
    if store:
        thread_logger.info("[Th %s]: Saving ...%s " % (chunk, func.__name__))
        majiq_io.dump_bin_file(results, "%s/%s_th%s.%s.pickle" % (tempdir, name, chunk, func.__name__))