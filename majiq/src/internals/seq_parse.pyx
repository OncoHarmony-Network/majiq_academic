from majiq.src.internals.io_bam cimport IOBam
from majiq.src.internals.grimoire cimport Junction, Gene, Exon, LSV, detect_lsvs, boostrap_samples
#
import cython
from majiq.src.constants import *
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.list cimport list as clist
from libcpp.vector cimport vector

import majiq.src.io as majiq_io
from cython.parallel import prange
from numpy.random  import choice
from scipy.stats import nbinom, poisson
import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_t
ctypedef vector[Gene *] gene_vect

# cdef int __mark_stacks(np.ndarray[np.float_t, ndim=2] junctions, float fitfunc_r, float pvalue_limit) :
#
#     cdef np.ndarray[np.float_t, ndim=2] pvalues
#     cdef np.ndarray[np.float_t, ndim=2] mean_rest
#     cdef np.ndarray[np.int_t, ndim=2] denom
#     cdef float r, p
#
#     if pvalue_limit <= 0:
#         return 0
#     denom = np.count_nonzero(junctions, axis=1)[:, None] - (junctions > 0)
#     with np.errstate(divide='ignore',invalid='ignore'):
#         mean_rest = (junctions.sum(axis=1)[:, None] - junctions) / denom
#         mean_rest[np.isnan(mean_rest)] = 0.5
#     if fitfunc_r >0:
#         r = 1/fitfunc_r
#         p = r/(mean_rest + r)
#         pvalues = 1 - nbinom.cdf(junctions, r, p)
#     else:
#         pvalues = 1 - poisson.cdf(junctions, mean_rest)
#     junctions[pvalues<pvalue_limit] = 0
#
#
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# cdef np.ndarray[DTYPE_t, ndim=2] _bootstrap_samples(np.ndarray[DTYPE_t, ndim=2] junction_list, int m, int k):
#     """Given the filtered reads, bootstrap samples from every junction
#     :param junction_list:
#     :param m:
#     :param k:
#     :param discardzeros:
#     :param trimborder:
#     :param fitted_one_over_r:
#     :return:
#
#     """
#
#     cdef np.ndarray[DTYPE_t, ndim=2] all_samples = np.zeros(shape=(junction_list.shape[0], m), dtype=np.float)
#     cdef int npos_mult
#     cdef int iternumber
#     cdef float r = 0
#     cdef np.ndarray[DTYPE_t, ndim=1] junction, km_samples_means
#     cdef int i
#
#     for i in range(junction_list.shape[0]):
#         junction = junction_list[i][junction_list[i] > 0]
#         npos_mult = np.count_nonzero(junction)
#         if npos_mult > 0:
#
#             all_samples[i, :m]  = np.reshape(choice(junction, k*m), (m, k)).mean(axis=1) * npos_mult
#             # all_samples[i, m] = junction_list[i].sum()
#             # all_samples[i, m+1] = npos_mult
#
#     return all_samples


cdef int _read_junction(list row, map[string, Junction*] jjs, map[string, Exon*] exs,
                        int nsamples, unsigned int eff_len) except -1:
    cdef string key = ('%s-%s' % (row[0], row[1])).encode('utf-8')
    jjs[key] = new Junction(row[0], row[1], nsamples, eff_len) #, annot=bool(row[2]))


cdef int _read_exon(list row, map[string, Junction*] jjs, map[string, Exon*] exs) except -1:
    cdef string key = ('%s-%s' % (row[0], row[1])).encode('utf-8')
    exs[key] = new Exon(row[0], row[1], bool(row[2]))


cdef int _pass_ir(list row, str gne_id, map[string, Junction*] jjs, map[string, Exon*] exs) except -1:
    pass



cdef from_matrix_to_objects(Gene * gg, object elements, int nsamples, unsigned int eff_len):

    cdef dict func_list
    cdef list elem

  #  func_list = {EX_TYPE: _read_exon, IR_TYPE: _pass_ir, J_TYPE: _read_junction}

    for elem in elements:
        pass
        if elem[3] == EX_TYPE:
            _read_exon(elem, gg.junc_map, gg.exon_map)
        elif elem[3] == J_TYPE:
            _read_junction(elem, gg.junc_map, gg.exon_map, nsamples, eff_len)



cdef int _gene_analysis(vector[pair[string, string]] list_pair_files, map[string, int] strandness, int nsamples,
                        gene_vect gene_list, int ksamples, int msamples,
                        unsigned int min_experiments, unsigned int eff_len, int minpos, int minreads) nogil:

    cdef int j, count=0
    cdef clist[LSV*] out_lsvlist
    cdef IOBam c_iobam
    cdef float[:, :] boots
    cdef float * boots_ptr;
    cdef LSV * lsvObj
    cdef Gene* gg

    for j in range(nsamples):
        c_iobam = IOBam(list_pair_files[j].first, strandness[list_pair_files[j].first], eff_len, nsamples, j)
        # for gg in gene_list:
        #     # gg.print_gene()
        c_iobam.find_junctions_from_region(gene_list)


    for gg in gene_list:
        gg.detect_exons()
        #TODO: IR detection for later
        count = count + detect_lsvs(out_lsvlist, gg, min_experiments, eff_len, minpos, minreads)
        for j in range(nsamples):
            for lsvObj in out_lsvlist:
                boots_ptr = boostrap_samples(lsvObj, msamples, ksamples, j, eff_len)
        del gg
        # with gil: print ('KK4')

    gene_list.clear()
    return 0

cdef _extract_junctions(list file_list, object genes_dict, object elem_dict, conf, logger):

    cdef int n = len(genes_dict)
    cdef int nbatches
    cdef int nthreads = min(conf.nthreads, len(conf.sam_list))
    cdef int nsamples = len(file_list)
    cdef int minpos = conf.minpos
    cdef int minreads = conf.minreads
    cdef int i, j
    cdef int k=conf.k, m=conf.m
    cdef float pvalue_limit=conf.pvalue_limit
    cdef unsigned int min_experiments = 1 if conf.min_exp == -1 else conf.min_exp
    cdef unsigned int eff_len = conf.readLen - 2*MIN_BP_OVERLAP
    cdef Gene * gg

    cdef map[string, int] strandness
    cdef vector[gene_vect] gene_list
    cdef vector[pair[string, string]] list_pair_files
    cdef batchsize = 200;
    cdef char st = '+'
    cdef int bidx = -1
    cdef gene_vect vt
    nbatches = np.ceil(n/batchsize)

    for i, (gne_id, gene_obj) in enumerate(genes_dict.items()):
        gg = new Gene(gne_id.encode('utf-8'), gene_obj['name'].encode('utf-8'), gene_obj['chromosome'].encode('utf-8'),
                      st, gene_obj['start'], gene_obj['end'])
        from_matrix_to_objects(gg, elem_dict[gne_id], nsamples, eff_len)
        # gg.print_gene()

                      # gene_obj['strand'], gene_obj['start'], gene_obj['end'])
        if i%batchsize == 0 :
            j = 100 if i+100 < n else n-i
            bidx += 1
            vt = gene_vect()
            gene_list.push_back(vt)
            # print("PPP")
        # print("KKK", gne_id, "MM", bidx,"LLL", i%batchsize, "+++", nbatches)
        gene_list[bidx].push_back(gg)

    for exp_name, is_junc_file, name in file_list:
        cs1 = ('%s/%s.%s' % (conf.sam_dir, exp_name, SEQ_FILE_FORMAT)).encode('utf-8')
        list_pair_files.push_back((pair[string, string])(cs1, name.encode('utf-8')))
        strandness[cs1] = conf.strand_specific[exp_name]

    for i in prange(nbatches, nogil=True, num_threads=nthreads):
        with gil:
            logger.info("%s/%s" %(i, nbatches))
        _gene_analysis(list_pair_files, strandness, nsamples, gene_list[i], k, m,
                           min_experiments, eff_len, minpos, minreads)

## OPEN API FOR PYTHON

def find_new_junctions(list file_list, int chunk, object slf, object conf, object logger):
    _extract_junctions(file_list,  slf.genes_dict, slf.elem_dict, conf, logger)

