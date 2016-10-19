import argparse
import cPickle
import glob
import multiprocessing
import os
from multiprocessing import Pool
from multiprocessing.queues import JoinableQueue

import h5py

from voila.utils import utils_voila

logger = utils_voila.get_logger('convert.log')


def convert(pickle_files):
    def extension(path):
        return os.path.splitext(path)[1]

    def load_pickle(filename):
        logger.info('loading {0}...'.format(filename))
        with open(filename, 'r') as p:
            return cPickle.load(p)

    def splice_graph(pickle_file):
        pickle = load_pickle(pickle_file)
        logger.info('converting splice graph {0}'.format(pickle_file))
        with h5py.File(pickle_file + '.h5', 'w') as h:
            for gs in pickle:
                gs.to_hdf5(h)

    def pickle(pickle_file):
        pickle = load_pickle(pickle_file)
        logger.info('converting pickle {0}'.format(pickle_file))
        with h5py.File(pickle_file + '.h5', 'w') as h:
            pickle.to_hdf5(h)

    def worker():
        while True:
            pickle_file = queue.get()
            ext = extension(pickle_file)
            convert_func[ext](pickle_file)
            logger.info('done {0}.'.format(pickle_file))
            queue.task_done()

    convert_func = {'.splicegraph': splice_graph, '.pickle': pickle}
    files_count = 0
    queue = JoinableQueue()

    # producer
    for pickle_file in pickle_files:
        pickle_file_norm = os.path.abspath(os.path.expanduser(pickle_file))
        pickle_glob = (f for f in glob.glob(pickle_file_norm) if extension(f) in convert_func)
        for f in pickle_glob:
            files_count += 1
            queue.put(f)

    process_count = min(files_count, multiprocessing.cpu_count())
    pool = Pool(process_count, worker, maxtasksperchild=1)

    pool.close()
    queue.join()
    queue.close()


def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert pickle and splicegraph files to their HDF5 counterparts.')
    parser.add_argument('pickle_files', nargs='+', help='List of files to convert.  Globbing is supported.')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    convert(args.pickle_files)
