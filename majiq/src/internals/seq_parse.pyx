from majiq.src.internals.io_bam cimport IOBam, gene_vect
# from majiq.src.internals.interval cimport Interval, ITNode, insert
from majiq.src.internals.grimoire cimport Junction, Gene, Exon, LSV, detect_lsvs, boostrap_samples, sortGeneList
#
import cython
from majiq.src.constants import *
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.list cimport list as clist
from libcpp.vector cimport vector

cimport majiq.src.io as majiq_io
from cython.parallel import prange
from libc.stdlib cimport malloc, free
from numpy.random  import choice
from scipy.stats import nbinom, poisson
import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_t




cdef _find_junctions(list file_list, vector[Gene*] gene_vec,  object conf, object logger):

    cdef int n = gene_vec.size()
    cdef int nthreads = conf.nthreads
    cdef int nsamples = len(file_list)
    cdef int minpos = conf.minpos
    cdef int minreads = conf.minreads
    cdef int i, j
    cdef int k=conf.k, m=conf.m
    cdef float pvalue_limit=conf.pvalue_limit
    cdef unsigned int min_experiments = 1 if conf.min_exp == -1 else conf.min_exp
    cdef unsigned int eff_len = conf.readLen - 2*MIN_BP_OVERLAP

    cdef map[string, gene_vect] gene_list
    cdef Gene * gg


    # cdef char st = '+'
    cdef clist[LSV*] out_lsvlist
    cdef IOBam c_iobam
    # cdef np.ndarray buff


    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
    cdef np.ndarray junc_ids
    # cdef Interval intv
    cdef vector[Junction *] jvec
    # cdef int counter = 0

    cdef map[string, int] strandness
    cdef vector[pair[string, string]] list_pair_files
    cdef map[string, unsigned int] j_ids

    for exp_name, is_junc_file, name in file_list:
        cs1 = ('%s/%s.%s' % (conf.sam_dir, exp_name, SEQ_FILE_FORMAT)).encode('utf-8')
        list_pair_files.push_back((pair[string, string])(cs1, name.encode('utf-8')))
        strandness[cs1] = conf.strand_specific[exp_name]

    sortGeneList(gene_vec)

    for gobj in gene_vec:
        if gene_list.count(gobj.chromosome) == 0:
            gene_list[gobj.chromosome] = gene_vect()
        gene_list[gobj.chromosome].push_back(gobj)


    for j in range(nsamples):

        logger.info('Reading file %s' %(file_list[j][0]))
        with nogil:
            c_iobam = IOBam(list_pair_files[j].first, strandness[list_pair_files[j].first], eff_len,
                            nthreads, gene_list)
            c_iobam.ParseJunctionsFromFile()
            njunc = c_iobam.get_njuncs()

        boots = np.zeros(shape=(njunc, m), dtype=np.float32)
        junc_ids = np.chararray(shape=njunc, itemsize=250, order='C')

        with nogil:
            c_iobam.boostrap_samples(m, k, <np.float32_t *> boots.data)
        j_ids = c_iobam.get_junc_map()
        jvec = c_iobam.get_junc_vec()
        for it in j_ids:
            junc_ids[it.second] = it.first
            jvec[it.second].clear_nreads()

        logger.info('Done Reading file %s' %(file_list[j][0]))
        vals = {'bootstrap': boots, 'junc_ids': junc_ids}

        out_file = "%s/%s.juncs" % (conf.outDir, file_list[j][0])
        with open(out_file, 'w+b') as ofp:
            np.savez(ofp, **vals)


    for i in prange(n, nogil=True, num_threads=nthreads):
        gg = gene_vec[i]
        gg.detect_exons()
        #TODO: IR detection for later
        # detect_lsvs(out_lsvlist, gg, min_experiments, minpos, minreads)


    #Read bootstrap and fill the



        # vals.clear()

## OPEN API FOR PYTHON
cdef _core_build(str transcripts, list file_list, object conf, object logger):
    cdef vector[Gene *] gene_vec
    majiq_io._read_gff(transcripts, gene_vec, logger)
    _find_junctions(file_list, gene_vec, conf, logger)

cpdef core_build(str transcripts, list file_list, object conf, object logger):
    _core_build(transcripts, file_list, conf, logger)





