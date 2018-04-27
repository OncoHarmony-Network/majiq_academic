from majiq.src.internals.io_bam cimport IOBam, gene_vect
# from majiq.src.internals.interval cimport Interval, ITNode, insert
from majiq.src.internals.grimoire cimport Junction, Gene, Exon, LSV, Jinfo, detect_lsvs, boostrap_samples, sortGeneList
#
import cython
from majiq.src.constants import *
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.list cimport list as clist
from libcpp.vector cimport vector
cimport openmp
cimport majiq.src.io as majiq_io
from cython.parallel import prange

#from cython cimport parallel
from libc.stdlib cimport calloc, free
from numpy.random  import choice
from scipy.stats import nbinom, poisson
import numpy as np
cimport numpy as np



cdef _store_junc_file(np.ndarray boots, list junc_ids, str experiment_name, str outDir):

    cdef str out_file = "%s/%s.juncs" % (outDir, experiment_name)
    cdef dict vals = {'bootstrap': boots}
    dt = np.dtype('S250, u4, u4, f4, f4')
    vals['junc_info'] = np.array(junc_ids, dtype=dt)
    with open(out_file, 'w+b') as ofp:
        np.savez(ofp, **vals)

from cpython cimport PyObject, Py_INCREF


cdef int _output_lsv_file_single(vector[LSV*] out_lsvlist, str experiment_name, str outDir, int nthreads,
                                 unsigned int msamples):
    cdef dict cov_dict = {}
    cdef int nlsv = out_lsvlist.size()
    cdef str out_file, junc_file
    cdef list type_list = []
    cdef LSV* lsv_ptr
    cdef int njunc
    cdef Jinfo jobj_ptr
    cdef map[string, Jinfo] tlb_juncs
    cdef int i, j, junc_idx
    cdef list junc_info = []
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
    cdef np.ndarray junc_ids
    cdef object all_juncs
    cdef string key
    cdef string lsvid


    junc_file = "%s/%s.juncs" % (outDir, experiment_name)
    out_file = "%s/%s.majiq" % (outDir, experiment_name)
    with open(junc_file, 'rb') as fp:
        all_juncs = np.load(fp)
        boots = all_juncs['bootstrap']
        junc_ids = all_juncs['junc_info']

        ''' If we move towards this solution we can remove elements in Jinfo and make it smaller'''
        for i in range(junc_ids.shape[0]):
            tlb_juncs[junc_ids[i][0]] = Jinfo(i, junc_ids[i][1], junc_ids[i][2], junc_ids[i][3], junc_ids[i][4])
        del junc_ids

        for j in prange(nlsv, nogil=True, num_threads=nthreads):
            lsv_ptr = out_lsvlist[j]
            njunc = lsv_ptr.get_num_junctions()
            lsvid = lsv_ptr.get_id()
            with gil:
                cov_dict[lsvid.decode('utf-8')] = np.zeros(shape=(njunc, msamples), dtype=np.float32)
                type_list.append((lsvid.decode('utf-8'), lsv_ptr.get_type()))
            junc_idx = 0
            for junc in lsv_ptr.get_junctions():
                key = junc.get_key(lsv_ptr.get_gene())
                if tlb_juncs.count(key) > 0 :
                    jobj_ptr = tlb_juncs[key]
                    with gil:
                        cov_dict[lsvid.decode('utf-8')][junc_idx] = boots[tlb_juncs[key].index]
                        junc_info.append((lsvid.decode('utf-8'), junc.get_start(), junc.get_end(),
                                          jobj_ptr.sreads, jobj_ptr.npos))
                junc_idx = junc_idx + 1

        majiq_io.dump_lsv_coverage(out_file, cov_dict, type_list, junc_info, experiment_name)
        tlb_juncs.clear()

    return len(type_list)

# cdef bint _lsv_gather(LSV* lsv_ptr, np.ndarray[np.float32_t, ndim=2, mode="c"] source,
#                       np.ndarray[np.float32_t, ndim=2, mode="c"] target, map[string, Jinfo] tlb_juncs,
#                       clist[Jinfo*] junc_info_list, unsigned int msamples):
#     cdef bint pres
#     with nogil:
#         pres = lsv_ptr.gather_lsv_info(<np.float32_t *> source.data, <np.float32_t *> target.data, junc_info_list,
#                                        tlb_juncs, msamples)
#     return pres
#
# cdef int _output_lsv_file(vector[LSV*] out_lsvlist, str experiment_name, str outDir, int nthreads,
#                        unsigned int msamples):
#
#     cdef dict cov_dict = {}
#     cdef int nlsv = out_lsvlist.size()
#     cdef str out_file, junc_file
#     cdef list type_list = []
#
#     cdef clist[Jinfo*] junc_info_list
#
#     cdef bint pres
#     cdef LSV* lsv_ptr
#     cdef int njunc
#     cdef Jinfo* jobj_ptr
#     cdef map[string, Jinfo] tlb_juncs
#     cdef int i, j
#
#     cdef list junc_info = []
#
#     cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
#     cdef np.ndarray val
#     cdef np.ndarray junc_ids
#     cdef object all_juncs
#
#
#
#
#     junc_file = "%s/%s.juncs" % (outDir, experiment_name)
#     out_file = "%s/%s.majiq" % (outDir, experiment_name)
#     with open(junc_file, 'rb') as fp:
#         all_juncs = np.load(fp)
#         boots = all_juncs['bootstrap']
#         junc_ids = all_juncs['junc_info']
#         print ("##1", type(all_juncs['bootstrap']), np.isfortran(all_juncs['bootstrap']))
#         print ("##2", type(boots.data), np.isfortran(boots))
#
#         for i in range(junc_ids.shape[0]):
#             tlb_juncs[junc_ids[i][0]] = Jinfo(i, junc_ids[i][1], junc_ids[i][2], junc_ids[i][3], junc_ids[i][4])
#         del junc_ids
#         for j in prange(nlsv, nogil=True, num_threads=nthreads):
#
#             lsv_ptr = out_lsvlist[j]
#             njunc = lsv_ptr.get_num_junctions()
#             with gil:
#                 cov_dict[lsv_ptr.get_id().decode('utf-8')] = np.zeros(shape=(njunc, msamples), dtype=np.float32)
#                 pres =  _lsv_gather(out_lsvlist[j], boots, cov_dict[lsv_ptr.get_id().decode('utf-8')], tlb_juncs,
#                                     junc_info_list, msamples)
#             if not pres:
#                 continue
#             with gil:
#                 type_list.append((lsv_ptr.get_id(),lsv_ptr.get_type()))
#                 for jobj_ptr in junc_info_list:
#                     junc_info.append((lsv_ptr.get_id(), jobj_ptr.start, jobj_ptr.end, jobj_ptr.sreads, jobj_ptr.npos))
#
#         #     lsv_ptr = out_lsvlist[j]
#         #     njunc = lsv_ptr.get_num_junctions()
#         #     with gil:
#         #         cov_dict[lsv_ptr.get_id().decode('utf-8')] = np.zeros(shape=(njunc, msamples), dtype=np.float32)
#         #         val_data =
#         #         # val = np.zeros(shape=(njunc, msamples), dtype=np.float32)
#         #     pres = lsv_ptr.gather_lsv_info(<np.float32_t *> boots.data, <np.float32_t *> val.data, junc_info_list,
#         #                                    tlb_juncs, msamples)
#         #
#         #         # val_data = <np.float32_t *> pp
#         #         # val = np.zeros(shape=(njunc, msamples), dtype=np.float32)
#         #         # print(experiment_name, lsv_ptr.get_id(), njunc, val.shape, lsv_ptr.get_type())
#         #     # with gil:
#         #     #     print(pres, experiment_name, lsv_ptr.get_id(), njunc, val.shape, lsv_ptr.get_type())
#         #
#         #     if not pres:
#         #         continue
#         #     with gil:
#         #
#         #
#         #         # cov_dict[lsv_ptr.get_id().decode('utf-8')] = np.ndarray(shape=(njunc, msamples), dtype=np.float32,
#         #         #                                                         buffer=val_data, order='C' )
#         #         Py_INCREF(val.data)
#         #         cov_dict[lsv_ptr.get_id().decode('utf-8')] = val
#         #         print("PPPPP", id(cov_dict[lsv_ptr.get_id().decode('utf-8')].data), id(val.data))
#         #         type_list.append((lsv_ptr.get_id(),lsv_ptr.get_type()))
#         #         for jobj_ptr in junc_info_list:
#         #             junc_info.append((lsv_ptr.get_id(), jobj_ptr.start, jobj_ptr.end, jobj_ptr.sreads, jobj_ptr.npos))
#         #
#         majiq_io.dump_lsv_coverage(out_file, cov_dict, type_list, junc_info, experiment_name)
#         tlb_juncs.clear()
#
#     return len(type_list)


cdef _find_junctions(list file_list, vector[Gene*] gene_vec,  object conf, object logger):

    cdef int n = gene_vec.size()
    cdef int nthreads = conf.nthreads
    cdef int nsamples = len(file_list)
    cdef unsigned int minpos = conf.minpos
    cdef unsigned int minreads = conf.minreads
    cdef unsigned int denovo_thresh= conf.min_denovo
    cdef int i, j
    cdef int k=conf.k, m=conf.m
    cdef float pvalue_limit=conf.pvalue_limit
    cdef unsigned int min_experiments = 1 if conf.min_exp == -1 else conf.min_exp
    cdef unsigned int eff_len = conf.readLen - 2*MIN_BP_OVERLAP
    cdef int strandness, njunc
    cdef list group_list
    cdef int last_it_grp

    cdef Gene * gg
    cdef vector[LSV*] out_lsvlist
    cdef IOBam c_iobam
    cdef string name, bamfile
    cdef str tmp_str
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
    cdef list junc_ids
    cdef int nlsv

    cdef int* jvec


    cdef map[string, gene_vect] gene_list
    cdef map[string, unsigned int] j_ids
    cdef map[string, int*] junc_summary
    cdef pair[string, unsigned int] it

    # print("SIZE", gene_vec.size())
    sortGeneList(gene_vec)

    for gobj in gene_vec:
        if gene_list.count(gobj.get_chromosome()) == 0:
            gene_list[gobj.get_chromosome()] = gene_vect()
        gene_list[gobj.get_chromosome()].push_back(gobj)

    for tmp_str, group_list in conf.tissue_repl.items():
        name = tmp_str.encode('utf-8')
        last_it_grp = group_list[len(group_list) - 1]
        for j in group_list:
            logger.info('Reading file %s' %(file_list[j][0]))
            bamfile = ('%s/%s.%s' % (conf.sam_dir, file_list[j][0], SEQ_FILE_FORMAT)).encode('utf-8')
            strandness = conf.strand_specific[file_list[j][0]]

            with nogil:
                c_iobam = IOBam(bamfile, strandness, eff_len, nthreads, gene_list)
                c_iobam.ParseJunctionsFromFile()
                njunc = c_iobam.get_njuncs()

            boots = np.zeros(shape=(njunc, m), dtype=np.float32)
            with nogil:
                c_iobam.boostrap_samples(m, k, <np.float32_t *> boots.data)
                j_ids = c_iobam.get_junc_map()
                jvec = c_iobam.get_junc_vec_summary()

            logger.info("Update flags")
            for i in prange(n, nogil=True, num_threads=nthreads):
                gg = gene_vec[i]
                gg.update_junc_flags(eff_len, (j==last_it_grp), minreads, minpos, denovo_thresh, min_experiments)

            logger.info("Done Update flags")
            junc_ids = [0] * njunc
            for it in j_ids:
                # print(file_list[j][0], j, last_it_grp, (j==last_it_grp), it.second, it.first)
                tmp_str = it.first.decode('utf-8').split(':')[2]
                start, end = (int(xx) for xx in tmp_str.split('-'))
                junc_ids[it.second] = (it.first.decode('utf-8'), start, end, jvec[it.second], jvec[it.second + njunc])

            logger.info('Done Reading file %s' %(file_list[j][0]))
            _store_junc_file(boots, junc_ids, file_list[j][0], conf.outDir)

    logger.info("Detecting LSVs ngenes:%s " % n)
    for i in prange(n, nogil=True, num_threads=nthreads):
        gg = gene_vec[i]
        gg.detect_exons()
        # with gil:
        #     gene_to_splicegraph(gg, conf)
        # gg.print_exons()

        #TODO: IR detection for later
        nlsv = detect_lsvs(out_lsvlist, gg)
        with gil:
            print("GENE: %s %s" % (gg.get_id(), nlsv))

    logger.info("Generating Splicegraph")
    # for i in range(n):
    #     gene_to_splicegraph(gene_vec[i], conf)

    logger.info("%s LSV found" % out_lsvlist.size())
    # #Read bootstrap and fill the
    logger.info('Storing LSVs files')
    for j in range(nsamples):
    #
        cnt = _output_lsv_file_single(out_lsvlist, file_list[j][0], conf.outDir, nthreads, m)
        logger.info('%s: %d LSVs' %(file_list[j][0], cnt))


from voila.api import SpliceGraph
from majiq.src.config import Config

cdef init_splicegraph(filename):
    majiq_config = Config()

    # erase splice graph file
    with SpliceGraph(filename, delete=True) as sg:
        sg.experiment_names = majiq_config.exp_list
        sg.genome = majiq_config.genome

cdef gene_to_splicegraph(Gene * gne, majiq_config):

    cdef pair[string,Junction *] jj_pair
    cdef Junction * jj
    cdef pair[string,Exon *] ex_pair
    cdef Exon * ex

    alt_empty_starts = []
    alt_empty_ends = []
    gne_id = gne.get_id()
    # print('GENEID:', gne_id)
    with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir)) as sg:
        sg.gene(gne_id).add(name=gne.get_name(), strand=gne.get_strand(), chromosome=gne.get_chromosome())

        for jj_pair in gne.junc_map_:
            jj = jj_pair.second
            if jj.get_start() == FIRST_LAST_JUNC:
                alt_empty_starts.append(jj.get_end())
                continue

            if jj.get_end() == FIRST_LAST_JUNC:
                alt_empty_ends.append(jj.get_start())
                continue

            # TODO: add transcripts
            sg.junction(gne_id, jj.get_start(), jj.get_end()).add(annotated=jj.get_annot(),
                                                                  intron_retention=jj.get_intronic())

        for ex_pair in gne.exon_map_:
            ex = ex_pair.second
        #for ex in sorted(exon_dict, key=lambda x: (x.start, x.end)):
            if ex.intron_:
                continue

            covered = False
            alt_start = []
            for jj in ex.ib:
                covered = covered or (jj.get_bld_fltr())
                if jj.get_end() in alt_empty_starts:
                    alt_start.append(jj.get_end())

            alt_ends = []
            for jj in ex.ob:
                covered = covered or (jj.get_bld_fltr())
                if jj.get_start() in alt_empty_ends:
                    alt_ends.append(jj.get_start())

            extra_coords = []
            if ex.annot_:
                if ex.start_ < ex.db_start_:
                    extra_coords.append([ex.start_, ex.db_start_ - 1])
                if ex.end_ > ex.db_end_:
                    extra_coords.append([ex.db_end_ + 1, ex.end_])

            sg.exon(gne_id, ex.start_, ex.end_).add(coords_extra=extra_coords, intron_retention=False,
                                                    annotated=ex.annot_, alt_starts=alt_start, alt_ends=alt_ends)

        # for info in list_introns:
        #     if info.skip:
        #         continue
        #
        #     sg.exon(gne_id, info.start, info.end).add(annotated=info.annot, intron_retention=True)



## OPEN API FOR PYTHON
cdef _core_build(str transcripts, list file_list, object conf, object logger):
    cdef vector[Gene *] gene_vec
    majiq_io.read_gff(transcripts, gene_vec, logger)
    init_splicegraph(get_builder_splicegraph_filename(conf.outDir))
    _find_junctions(file_list, gene_vec, conf, logger)

cpdef core_build(str transcripts, list file_list, object conf, object logger):
    _core_build(transcripts, file_list, conf, logger)





