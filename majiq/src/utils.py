import logging
import os
import random
import resource
import sys

import h5py
import numpy as np

import majiq.grimoire.lsv as majiq_lsv
import majiq.src.config as majiq_config
import majiq.grimoire.junction as majiq_junction
import majiq.src.io_utils as majiq_io_utils




class Writer(object):
    """Create an object with a write method that writes to a
    specific place on the screen, defined at instantiation.

    This is the glue between blessings and progressbar.
    """
    def __init__(self, location):
        """
        Input: location - tuple of ints (x, y), the position
                        of the bar in the terminal
        """
        self.location = location

    def write(self, string):
        with majiq_config.term.location(*self.location):
            print(string)

def monitor(msg):
    print "MONITOR", msg, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000, 'MB'
    sys.stdout.flush()


def create_if_not_exists(my_dir, logger=False):
    """Create a directory path if it does not exist"""
    try:
        if logger:
            logger.debug("\nCreating directory %s..." % my_dir)
        os.makedirs(my_dir)
    except OSError:
        if logger:
            logger.debug("\nDirectory %s already exists..." % my_dir)


def get_logger(logger_name, silent=False, debug=False, child=False):
    """
    Returns a logger instance. verbose = False will silence the logger, debug will give 
    more information intended for debugging purposes.
    """
    logging_format = "%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s"
    logger = logging.getLogger(logger_name)
    formatter = logging.Formatter(logging_format)

    fileHandler = logging.FileHandler(logger_name, mode='w')
    fileHandler.setFormatter(formatter)

    streamHandler = logging.StreamHandler()
    streamHandler.setFormatter(formatter)

    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    logger.addHandler(fileHandler)
    logger.addHandler(streamHandler)

    return logger


def get_fitfunc_from_rnafile(path):
    with h5py.File(path) as p:
        res = p.attrs['one_over_r']
    return res


"""
Majiq file generation
"""


def send_output(lsv_list, non_as, temp_dir, out_queue, chnk, mlock):

## per bam file
    # for name, ind_list in majiq_config.tissue_repl.items():
    #
    #     njuncs = len(non_as[name])
    #     t_juncs = float(majiq_config.nrandom_junctions) / majiq_config.num_final_chunks
    #     prb = min(1.0, float(t_juncs) / njuncs) * 100
    #     kk = np.random.choice(100, njuncs)
    #     indx = np.arange(njuncs)[kk <= prb]
    #     r_junctions = np.array(list(non_as[name]))[indx]
    # #
    #     for idx, exp_idx in enumerate(ind_list):
    #
    #         for jix, jn in enumerate(r_junctions):
    #             out_queue.put([1, jn.get_coverage(exp_idx), exp_idx], block=True)

    out_queue.put([3, chnk, -1], block=True)
    mlock.acquire()
    mlock.release()
    # out_queue.close()


def prepare_lsv_table(lsv_list, non_as, temp_dir):
    """

    :param lsv_list:
    :param non_as:
    :param temp_dir:
    """
    for name, ind_list in majiq_config.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):
            fname = "%s/%s.majiq.hdf5" % (temp_dir, majiq_config.exp_list[exp_idx])
            f = h5py.File(fname, 'w', compression='gzip', compression_opts=9)
            as_table = f.create_group('LSVs')
            for iix, lsv in enumerate(lsv_list[name]):
                lsv.hdf5_lsv(as_table, exp_idx)

            non_as_table = f.create_group('const')
            for jix, jn in enumerate(non_as[name]):
                jn.hdf5_junction(non_as_table, exp_idx)
            f.close()


def merge_and_create_majiq_file(exp_idx, pref_file):
    """
    :param exp_idx: Index of experiment in config file
    :param pref_file: Prefix for the majiq name
    """

    experiment = majiq_config.exp_list[exp_idx]
    all_visual = []
    for chnk in range(majiq_config.num_final_chunks):
        temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
        temp_filename = '%s/%s.splicegraph.pkl' % (temp_dir, experiment)
        if os.path.exists(temp_filename):
            visual_gene_list = majiq_io_utils.load_bin_file(temp_filename)
            all_visual.append(visual_gene_list)
    fname = '%s/%s%s.splicegraph' % (majiq_config.outDir, pref_file, experiment)
    visual = np.concatenate(all_visual)
    majiq_io_utils.dump_bin_file(visual, fname)
    del all_visual
    del visual

    fname = '%s/%s%s.majiq' % (majiq_config.outDir, pref_file, experiment)
    of = h5py.File(fname, 'w', compression='gzip', compression_opts=9)

    of['experiment'] = experiment
    of['genome'] = majiq_config.genome
    of['num_reads'] = majiq_config.num_mapped_reads[exp_idx]
    as_table = of.create_group('LSVs')
    nonas_table = of.create_group('const')

    nat = []
    for chnk in range(majiq_config.num_final_chunks):
        temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
        filename = "%s/%s.majiq.hdf5" % (temp_dir, majiq_config.exp_list[exp_idx])
        if os.path.exists(filename):
            with h5py.File(filename) as temp_table:
                for kk in temp_table['LSVs'].keys():
                    h5py.h5o.copy(temp_table['LSVs'].id, kk, as_table.id, kk)
                    majiq_lsv.set_gc_factor(as_table[kk], exp_idx)
                for kk in temp_table['const'].keys():
                    nat.append([kk, temp_table.filename])

    clist = random.sample(nat, min(5000, len(nat)))
    for jnc in clist:
        with h5py.File(jnc[1]) as tt:
            h5py.h5o.copy(tt['const'].id, jnc[0], nonas_table.id, jnc[0])
            majiq_junction.set_gc_factor(nonas_table[jnc[0]], exp_idx)


def print_junc_matrices(mat, tlb=None, fp=None):
    """

    :param mat:
    :param tlb:
    :param fp:
    """
    if not fp is None:
        out = open('./junc_matrix.tab', 'a+')
    else:
        out = sys.stdout
    out.write("\n=== BEGIN %s === \n\n" % id)
    N, M = mat.shape
    header = [0] * N
    if not tlb is None:
        out.write("Nan\t")
        for ex, (p1, p2) in tlb.items():
            for n in p1:
                out.write("%d\t" % (ex + 1))
            for n in p2:
                header[n] = "%d" % (ex + 1)
    out.write("\n")
    for ii in np.arange(N):
        if not tlb is None:
            out.write("%s\t" % header[ii])
        for jj in np.arange(M):
            val = mat[ii, jj]
            out.write("%s\t" % val)
        out.write("\n")
    out.write("\n=== END %s === \n\n" % id)
    if fp is None:
        out.close()


def get_validated_pcr_lsv(candidates, out_dir):
    """

    :param candidates:
    :param out_dir:
    """
    pcr_list = []
    print "get_validated_pcr_lsv", len(candidates[0])
    for lsv in candidates[0]:
        if not lsv.has_pcr_score():
            continue
        alt_coord = lsv.exon.get_pcr_candidate()
        score = lsv.get_pcr_score()
        for jidx, jj in enumerate(lsv.junctions):
            if lsv.is_Ssource():
                excoord = jj.get_acceptor().get_coordinates()
            else:
                excoord = jj.get_donor().get_coordinates()
            if excoord[1] > alt_coord[0] and excoord[0] < alt_coord[1]:
                name = "%s#%s" % (lsv.id, jidx)
                pcr_lsv = [lsv.exon.get_pcr_name(), name, score]
                pcr_list.append(pcr_lsv)
                print "PCR", ' '.join(pcr_lsv)
    fname = '%s/pcr.pkl' % out_dir
    majiq_io_utils.dump_bin_file(pcr_list, fname)


# ANALYSIS FUNCTIONS

def analyze_denovo_junctions(genes, output):
    denovo_list = [[] for xx in range(majiq_config.num_experiments)]
    annot_list = [[] for xx in range(majiq_config.num_experiments)]

    for gg in genes:
        jlist = gg.get_all_junctions()
        for jj in jlist:
            for tissue, list_idx in majiq_config.tissue_repl.items():
                for exp_n in list_idx:
                    if jj.is_annotated():
                        annot_list[exp_n].append(jj)
                    else:
                        denovo_list[exp_n].append(jj)

    majiq_io_utils.dump_bin_file([majiq_config.tissue_repl, annot_list, denovo_list], output)


def histogram_for_exon_analysis(genes, output):
    denovo_list = []
    annotated_list = []

    for gg in genes:
        ex_list = gg.get_exon_list()
        for ex in ex_list:
            lngth = ex.get_length()
            if ex.is_annotated():
                annotated_list.append(lngth)
            else:
                denovo_list.append(lngth)

    majiq_io_utils.dump_bin_file([annotated_list, denovo_list], output)


def chunks2(l, n, extra):
    """Yield successive n-sized chunks from l.
    :param l: list to be split
    :param n: max length of chunks
    """
    try:
        idx = -1
        for i in range(0, len(l), n):
            idx += 1
            if extra is not None:
                yield (l[i:i+n], extra[idx])
            else:
                yield l[i:i+n]
    except:
        print "ERROR: extra value has incorrect size %s" %idx, extra
        raise


def chunks(l, n, extra):
    """Yield successive n-sized chunks from l.
    :param extra:
    :param l: list to be split
    :param n: max length of chunks
    """
    try:
        idx = -1
        rep_chunk = [0] * len(extra)
        for i in range(len(l)):
            idx += 1
            eidx = idx % len(extra)
            yield (l[i], extra[eidx], rep_chunk[eidx], n)
            rep_chunk[eidx] += 1

    except:
        print "ERROR: extra value has incorrect size %s" % idx, extra
        raise

import sys
from numbers import Number
from collections import Set, Mapping, deque

try: # Python 2
    zero_depth_bases = (basestring, Number, xrange, bytearray)
    iteritems = 'iteritems'
except NameError: # Python 3
    zero_depth_bases = (str, bytes, Number, range, bytearray)
    iteritems = 'items'

def getsize(obj_0):
    """Recursively iterate to sum size of object & members."""
    def inner(obj, _seen_ids = set()):
        obj_id = id(obj)
        if obj_id in _seen_ids:
            return 0
        _seen_ids.add(obj_id)
        size = sys.getsizeof(obj)
        if isinstance(obj, zero_depth_bases):
            pass # bypass remaining control flow and return
        elif isinstance(obj, (tuple, list, Set, deque)):
            size += sum(inner(i) for i in obj)
        elif isinstance(obj, Mapping) or hasattr(obj, iteritems):
            size += sum(inner(k) + inner(v) for k, v in getattr(obj, iteritems)())
        # Check for custom object instances - may subclass above too
        if hasattr(obj, '__dict__'):
            size += inner(vars(obj))
        if hasattr(obj, '__slots__'): # can have __slots__ with __dict__
            size += sum(inner(getattr(obj, s)) for s in obj.__slots__ if hasattr(obj, s))
        return size
    return inner(obj_0)