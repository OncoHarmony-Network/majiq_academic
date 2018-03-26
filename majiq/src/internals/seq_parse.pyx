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
from libc.stdlib cimport malloc, free
from numpy.random  import choice
from scipy.stats import nbinom, poisson
import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_t
ctypedef vector[Gene *] gene_vect


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

    # for j in range(nsamples):
    #     c_iobam = IOBam(list_pair_files[j].first, strandness[list_pair_files[j].first], eff_len, nsamples, j)
    #     # for gg in gene_list:
    #     #     # gg.print_gene()
    #     c_iobam.find_junctions_from_region(gene_list)


    for gg in gene_list:
        gg.detect_exons()
        #TODO: IR detection for later
        count = count + detect_lsvs(out_lsvlist, gg, min_experiments, eff_len, minpos, minreads)
        # for j in range(nsamples):
            # for lsvObj in out_lsvlist:
                # boots_ptr = boostrap_samples(lsvObj, msamples, ksamples, j, eff_len)
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

cdef _find_junctions(list file_list, object genes_dict, object elem_dict, conf, logger):

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
    # cdef Gene * gg

    cdef map[string, int] strandness
    cdef gene_vect gene_list
    cdef vector[pair[string, string]] list_pair_files
    cdef char st = '+'
    cdef clist[LSV*] out_lsvlist
    cdef IOBam c_iobam

    cdef np.ndarray buff
    cdef map[string, unsigned int] j_ids
    cdef np.npy_intp pn1[2] , pn2[2];

    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots ;
    cdef np.ndarray junc_ids ;



    # for i, (gne_id, gene_obj) in enumerate(genes_dict.items()):
    #     gg = new Gene(gne_id.encode('utf-8'), gene_obj['name'].encode('utf-8'), gene_obj['chromosome'].encode('utf-8'),
    #                   st, gene_obj['start'], gene_obj['end'])
    #     from_matrix_to_objects(gg, elem_dict[gne_id], nsamples, eff_len)
    #
    #     gene_list.push_back(gg)

    for exp_name, is_junc_file, name in file_list:
        cs1 = ('%s/%s.%s' % (conf.sam_dir, exp_name, SEQ_FILE_FORMAT)).encode('utf-8')
        list_pair_files.push_back((pair[string, string])(cs1, name.encode('utf-8')))
        strandness[cs1] = conf.strand_specific[exp_name]

    for j in range(nsamples):

        logger.info('Reading file %s' %(file_list[j][0]))
        with nogil:
            c_iobam = IOBam(list_pair_files[j].first, strandness[list_pair_files[j].first], eff_len, nthreads)
            c_iobam.ParseJunctionsFromFile()
            njunc = c_iobam.get_njuncs()

        boots = np.zeros(shape=(njunc, m), dtype=np.float32)
        junc_ids = np.chararray(shape=njunc, itemsize=250, order='C')

        with nogil:
            c_iobam.boostrap_samples(m, k, <np.float32_t *> boots.data)
        j_ids = c_iobam.get_junc_map()
        for it in j_ids:
            junc_ids[it.second] = it.first

        logger.info('Done Reading file %s' %(file_list[j][0]))
        vals = {'bootstrap': boots, 'junc_ids': junc_ids}

        out_file = "%s/%s.juncs" % (conf.outDir, file_list[j][0])
        with open(out_file, 'w+b') as ofp:
            np.savez(ofp, **vals)


        # vals.clear()

## OPEN API FOR PYTHON

def find_new_junctions2(list file_list, int chunk, object slf, object conf, object logger):
    _extract_junctions(file_list,  slf.genes_dict, slf.elem_dict, conf, logger)

def find_new_junctions(list file_list, int chunk, object slf, object conf, object logger):
    _find_junctions(file_list,  slf.genes_dict, slf.elem_dict, conf, logger)



