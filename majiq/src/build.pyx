import os
import sys
import psutil

from majiq.src.internals.grimoire cimport Junction, Gene, Exon, LSV, Jinfo, Intron
from majiq.src.internals.io_bam cimport IOBam, prepare_genelist, overGene_vect_t, free_genelist
from majiq.src.internals.grimoire cimport find_intron_retention, find_gene_from_junc, isNullJinfo, fill_junc_tlb
from majiq.src.internals.grimoire cimport key_format, free_JinfoVec, Gene_vect_t, free_lsvlist
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.polyfitnb cimport fit_nb
from majiq.src.config import Config
import majiq.src.logger as majiq_logger
cimport majiq.src.io as majiq_io
from majiq.src.constants import *
from voila.c.splice_graph_sql cimport open_db, close_db
from voila.c.splice_graph_sql cimport gene as sg_gene
from voila.c.splice_graph_sql cimport junction as sg_junction
from voila.c.splice_graph_sql cimport exon as sg_exon
from voila.c.splice_graph_sql cimport intron_retention as sg_intron_retention
from voila.c.splice_graph_sql cimport alt_start as sg_alt_start
from voila.c.splice_graph_sql cimport alt_end as sg_alt_end
from voila.c.splice_graph_sql cimport junction_reads as sg_junction_reads
from voila.c.splice_graph_sql cimport intron_retention_reads as sg_intron_retention_reads
from voila.api import SpliceGraph
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from cython.parallel import prange
import cython
from cython import parallel
import numpy as np
cimport numpy as np


cdef extern from "sqlite3.h":
    struct sqlite3

cdef int C_FIRST_LAST_JUNC = FIRST_LAST_JUNC


cdef _store_junc_file(np.ndarray boots, np.ndarray ir_cov, list junc_ids, str experiment_name, np.ndarray meta, str outDir):

    cdef str out_file = "%s/%s.%s" % (outDir, experiment_name, JUNC_FILE_FORMAT)
    cdef dict vals = {'bootstrap': boots, 'ir_cov': ir_cov}
    dt = np.dtype('S250, u4, u4, f4, f4, u4')
    vals['junc_info'] = np.array(junc_ids, dtype=dt)
    vals['meta'] = meta
    with open(out_file, 'w+b') as ofp:
        np.savez(ofp, **vals)


cdef void update_splicegraph_junction(sqlite3 *db, string gene_id, int start, int end, int nreads, string exp) nogil:
    if C_FIRST_LAST_JUNC != start and C_FIRST_LAST_JUNC != end:
        sg_junction_reads(db, nreads, exp, gene_id, start, end)

ctypedef vector[Jinfo *] jinfoptr_vec_t

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef int _output_majiq_file(vector[LSV*] lsvlist, map[string, overGene_vect_t] gList, map[string, int] j_tlb,
                            tuple fname, string outDir, sqlite3* db, unsigned int msamples,
                            bint irb, bint simpl, object logger, int nthreads) except -1:

    cdef unsigned int irbool, coord1, coord2, sreads, npos
    cdef unsigned int nlsv = lsvlist.size()
    cdef list cov_l, type_list = []
    cdef list junc_info = []
    cdef str out_file, junc_file
    cdef int njunc = 0
    cdef np.float32_t[:, :] boots
    # cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
    cdef np.ndarray junc_ids
    cdef int i, j, njlsv = j_tlb.size()
    cdef unsigned int junc_idx, m
    cdef jinfoptr_vec_t jobj_vec
    cdef Gene_vect_t gene_l
    cdef string key, chrom, lsvid, gid, jid
    cdef Jinfo* jobj_ptr
    cdef vector[np.float32_t] x
    cdef vector[vector[int]] tmp_juncinfo
    cdef vector[vector[np.float32_t]] tmp_boots
    cdef vector[Junction*] tmp_juncvec
    cdef char strand
    cdef string experiment_name = fname[0].encode('utf-8')

    cdef int thread_id = -1

    logger.info('DUMP file %s' % experiment_name)
    jobj_vec = jinfoptr_vec_t(njlsv)

    # sg_filename = get_builder_splicegraph_filename(outDir.decode('utf-8')).encode('utf-8')
    if fname[2]:
        junc_file = fname[1]
    else:
        junc_file = "%s/%s.%s" % (outDir.decode('utf-8'), experiment_name.decode('utf-8'), JUNC_FILE_FORMAT)
    out_file = "%s/%s.%s" % (outDir.decode('utf-8'), experiment_name.decode('utf-8'), MAJIQ_FILE_FORMAT)
    with open(junc_file, 'rb') as fp:
        junc_ids = np.load(fp)['junc_info']
    njunc = junc_ids.shape[0]

#TODO: CHECK THE  with gil statement if it is necessary
    for i in prange(njunc, nogil=True, num_threads=nthreads):
        gene_l = Gene_vect_t()
        with gil:
            jid     = junc_ids[i][0]
            coord1  = junc_ids[i][1]
            coord2  = junc_ids[i][2]
            sreads  = junc_ids[i][3]
            npos    = junc_ids[i][4]
            irbool  = junc_ids[i][5]
            chrom   = jid.split(b':')[0]
            strand  = <char> jid.split(b':')[1][0]

        find_gene_from_junc(gList, chrom, strand, coord1, coord2, gene_l, irbool, simpl)
        if irbool == 0:
            for gneObj in gene_l:
                update_splicegraph_junction(db, gneObj.get_id(), coord1, coord2, sreads, experiment_name)
                with gil:
                    key = key_format(gneObj.get_id(), coord1, coord2, False)
                    if j_tlb.count(key) > 0:
                        jobj_ptr = new Jinfo(i, sreads, npos)
                        jobj_vec[j_tlb[key]] = jobj_ptr

        elif irb:
            with gil:
                gid = b':'.join(jid.split(b':')[3:])
            for gneObj in gene_l:
                if gneObj.get_id() != gid:
                    continue
                irv = find_intron_retention(gneObj, coord1, coord2)
                for ir_ptr in irv:
                    sg_intron_retention_reads(db, sreads, experiment_name,  gneObj.get_id(),
                                              ir_ptr.get_start(), ir_ptr.get_end())
                    with gil:
                        key = key_format(gneObj.get_id(), ir_ptr.get_start(), ir_ptr.get_end(), True)
                        if j_tlb.count(key) > 0:
                            jobj_ptr = new Jinfo(i, sreads, npos)
                            jobj_vec[j_tlb[key]] = jobj_ptr

        gene_l.clear()

    del junc_ids
    logger.info("Create majiq file")

    with open(junc_file, 'rb') as fp:
        boots = np.load(fp)['bootstrap']
    cov_l = list()
    junc_info = []
    type_list = []

    for j in prange(nlsv, nogil=True, num_threads=nthreads):
        thread_id = parallel.threadid()
        lsv_ptr = lsvlist[j]
        njunc = lsv_ptr.get_num_variations()
        if njunc<2: continue
        lsvid = lsv_ptr.get_id()
        with gil:
            type_list.append((lsvid.decode('utf-8'), lsv_ptr.get_type()))

        tmp_juncinfo = vector[vector[int]](njunc, vector[int](4))
        tmp_boots    = vector[vector[np.float32_t]](njunc, vector[np.float32_t](msamples))
        tmp_juncvec  = lsv_ptr.get_junctions()

        for junc_idx in range(tmp_juncvec.size()):
            junc = tmp_juncvec[junc_idx]
            key = junc.get_key(lsv_ptr.get_gene())
            tmp_juncinfo[junc_idx][0] = junc.get_start()
            tmp_juncinfo[junc_idx][1] = junc.get_end()

            if j_tlb.count(key) > 0 and not isNullJinfo(jobj_vec[j_tlb[key]]):
                jobj_ptr = jobj_vec[j_tlb[key]]
                tmp_juncinfo[junc_idx][2] = jobj_ptr.sreads
                tmp_juncinfo[junc_idx][3] = jobj_ptr.npos
                for m in range(msamples):
                    tmp_boots[junc_idx][m] = boots[jobj_ptr.index][m]

        ir_ptr = lsv_ptr.get_intron()

        if irb and ir_ptr != <Intron * > 0:
            key = key_format(lsv_ptr.get_gene().get_id(), ir_ptr.get_start(), ir_ptr.get_end(), True)
            tmp_juncinfo[njunc-1][0] = ir_ptr.get_start()
            tmp_juncinfo[njunc-1][1] = ir_ptr.get_end()

            if j_tlb.count(key) > 0 and not isNullJinfo(jobj_vec[j_tlb[key]]):
                jobj_ptr = jobj_vec[j_tlb[key]]
                tmp_juncinfo[njunc-1][2] = jobj_ptr.sreads
                tmp_juncinfo[njunc-1][3] = jobj_ptr.npos
                for m in range(msamples):
                    tmp_boots[njunc-1][m] = boots[jobj_ptr.index][m]

        with gil:
            for i in range(njunc):
                junc_info.append((lsvid.decode('utf-8'), tmp_juncinfo[i][0], tmp_juncinfo[i][1],
                                  tmp_juncinfo[i][2],tmp_juncinfo[i][3]))
                tc = []
                for m in range(msamples):
                    tc.append(tmp_boots[i][m])
                cov_l.append(tuple(tc))

    logger.info("Dump majiq file")
    majiq_io.dump_lsv_coverage_mat(out_file, cov_l, type_list, junc_info, experiment_name.decode('utf-8'))
    free_JinfoVec(jobj_vec)
    nlsv = len(type_list)

    return nlsv


cdef _parse_junction_file(tuple filetp, map[string, Gene*]& gene_map, vector[string] gid_vec,
                          map[string, overGene_vect_t] gene_list, int min_experiments, bint reset, object conf,
                          object logger):

    cdef int nthreads = conf.nthreads
    cdef int strandness = conf.strand_specific[filetp[0]]
    cdef unsigned int minreads = conf.minreads
    cdef unsigned int minpos   = conf.minpos
    cdef unsigned int denovo_thresh= conf.min_denovo
    cdef bint denovo = conf.denovo
    cdef IOBam c_iobam
    cdef int njunc, i, j
    cdef np.ndarray junc_ids
    cdef object fp
    cdef Gene_vect_t gene_l
    cdef string key, chrom, lsvid, gid, jid
    # cdef int coord1, coord2, sreads, npos, strand
    cdef unsigned int irbool, coord1, coord2, sreads, npos
    cdef int n = gene_map.size()
    cdef char strand
    cdef bint bsimpl = (conf.simpl_psi >= 0)
    cdef bint ir = conf.ir
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] ir_cov
    cdef vector[np.float32_t] ir_vec
    cdef unsigned int eff_len = conf.readLen - 2*MIN_BP_OVERLAP + 1
    cdef np.float32_t min_ir_cov = conf.min_intronic_cov
    cdef np.float32_t ir_numbins = conf.irnbins
    cdef int jlimit

    c_iobam = IOBam(filetp[1].encode('utf-8'), strandness, eff_len, nthreads, gene_list, bsimpl)

    with np.load(filetp[1]) as fp:
        junc_ids = fp['junc_info']
        if ir:
            ir_cov = fp['ir_cov']
        jlimit = fp['meta'][0][2]
    njunc = junc_ids.shape[0]

    for j in prange(njunc, nogil=True, num_threads=nthreads):
        gene_l = Gene_vect_t()
        with gil:
            jid     = junc_ids[j][0]
            coord1  = junc_ids[j][1]
            coord2  = junc_ids[j][2]
            sreads  = junc_ids[j][3]
            npos    = junc_ids[j][4]
            irbool  = junc_ids[j][5]
            chrom   = jid.split(b':')[0]
            strand  = <char> jid.split(b':')[1][0]
            gid     = b'.'
            if irbool == 1 and not ir :
                continue
            elif ir and irbool == 1:
                gid = b':'.join(jid.split(b':')[3:])
                ir_vec = vector[np.float32_t](eff_len)
                for i in range(eff_len):
                    ir_vec[i] = ir_cov[j - jlimit][i]

                # logger.info("IR VEC: %s %s" %(eff_len, ir_vec.size()))
        c_iobam.parseJuncEntry(gene_list, gid, chrom, strand, coord1, coord2, sreads, minreads, npos, minpos,
                               denovo_thresh, denovo, gene_l, irbool==1, ir_vec, min_ir_cov, ir_numbins,
                               min_experiments, reset)


    # for i in prange(n, nogil=True, num_threads=nthreads):
    #     gg = gene_map[gid_vec[i]]
    #     gg.update_junc_flags(1, reset, minreads, 0, denovo_thresh, min_experiments, denovo)

    c_iobam.free_iobam()
    logger.info('Done Reading file %s' %(filetp[0]))


cdef _find_junctions(list file_list, map[string, Gene*]& gene_map, vector[string] gid_vec,
                     map[string, overGene_vect_t] gene_list, object conf, object logger):

    cdef int n = gene_map.size()
    cdef int nthreads = conf.nthreads

    cdef unsigned int minpos = conf.minpos
    cdef unsigned int minreads = conf.minreads
    cdef unsigned int denovo_thresh= conf.min_denovo
    cdef bint denovo = conf.denovo
    cdef np.float32_t min_ir_cov = conf.min_intronic_cov
    cdef int k=conf.k, m=conf.m
    cdef np.float32_t pvalue_limit=conf.pvalue_limit
    cdef unsigned int min_experiments
    cdef unsigned int eff_len = conf.readLen - 2*MIN_BP_OVERLAP + 1
    cdef bint ir = conf.ir
    cdef bint bsimpl = (conf.simpl_psi >= 0)
    cdef np.float32_t ir_numbins=conf.irnbins

    cdef int i, j
    cdef int strandness, njunc
    cdef list group_list
    cdef int last_it_grp

    cdef Gene * gg
    cdef IOBam c_iobam
    cdef string name, bamfile
    cdef str tmp_str
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] ir_raw_cov
    cdef list junc_ids
    cdef np.float32_t fitfunc_r
    cdef unsigned int jlimit
    cdef int* jvec
    cdef map[string, unsigned int] j_ids
    cdef pair[string, unsigned int] it

    for tmp_str, group_list in conf.tissue_repl.items():
        name = tmp_str.encode('utf-8')
        last_it_grp = group_list[len(group_list) - 1]
        min_experiments = conf.min_experiments[tmp_str]
        logger.info('Group %s, number of experiments: %s, minexperiments: %s' % (tmp_str,
                                                                                  len(group_list), min_experiments))
        for j in group_list:
            if file_list[j][2]:
                logger.info('Reading %s file %s' %(JUNC_FILE_FORMAT, file_list[j][1]))
                _parse_junction_file(file_list[j], gene_map, gid_vec,gene_list, min_experiments, (j==last_it_grp),
                                     conf, logger)
            else:
                logger.info('Reading %s file %s' %(SEQ_FILE_FORMAT, file_list[j][1]))
                bamfile = ('%s' % (file_list[j][1])).encode('utf-8')
                strandness = conf.strand_specific[file_list[j][0]]

                with nogil:
                    c_iobam = IOBam(bamfile, strandness, eff_len, nthreads, gene_list, bsimpl)
                    c_iobam.ParseJunctionsFromFile(False)
                    n_junctions = c_iobam.get_njuncs()
                    if ir:
                        with gil:
                            logger.info('Detect Intron retention %s' %(file_list[j][0]))
                        c_iobam.detect_introns(min_ir_cov, min_experiments, ir_numbins, (j==last_it_grp))
                    njunc = c_iobam.get_njuncs()
                    with gil:
                        logger.debug('Total Junctions and introns %s' %(njunc))

                if n_junctions == 0 or pvalue_limit <= 0:
                    if n_junctions == 0:
                        logger.warning('No junctions were found on sample %s' % bamfile)
                    fitfunc_r = 0
                else:
                    fitfunc_r = fit_nb(c_iobam.junc_vec, n_junctions, eff_len, nbdisp=0.1, logger=logger)

                boots = np.zeros(shape=(njunc, m), dtype=np.float32)
                with nogil:
                    c_iobam.boostrap_samples(m, k, <np.float32_t *> boots.data, fitfunc_r, pvalue_limit)
                    j_ids  = c_iobam.get_junc_map()
                    jvec   = c_iobam.get_junc_vec_summary()
                    jlimit = c_iobam.get_junc_limit_index()

                ir_raw_cov = np.zeros(shape=(njunc - jlimit, eff_len), dtype=np.float32)
                if ir:
                    with nogil:
                        c_iobam.get_intron_raw_cov(<np.float32_t *> ir_raw_cov.data)

                logger.debug("Update flags")
                for i in prange(n, nogil=True, num_threads=nthreads):
                    gg = gene_map[gid_vec[i]]
                    gg.update_junc_flags(eff_len, (j==last_it_grp), minreads, minpos, denovo_thresh, min_experiments, denovo)

                logger.debug("Done Update flags")
                junc_ids = [0] * njunc
                for it in j_ids:
                    tmp_str = it.first.decode('utf-8').split(':')[2]
                    start, end = (int(xx) for xx in tmp_str.split('-'))
                    junc_ids[it.second] = (it.first.decode('utf-8'), start, end, jvec[it.second], jvec[it.second + njunc],
                                           int(it.second>= jlimit))

                logger.info('Done Reading file %s' %(file_list[j][0]))
                dt = np.dtype('|S250, |S25, u4')
                meta = np.array([(file_list[j][0], VERSION, jlimit)], dtype=dt)
                _store_junc_file(boots, ir_raw_cov, junc_ids, file_list[j][0], meta, conf.outDir)
                c_iobam.free_iobam()
                del boots
                del ir_raw_cov


cdef init_splicegraph(string filename, object conf):

    with SpliceGraph(filename.decode('utf-8'), delete=True) as sg:
        sg.experiment_names = conf.exp_list
        sg.genome = conf.genome

cdef void gene_to_splicegraph(Gene * gne, sqlite3 * db) nogil:

    cdef pair[string, Junction *] jj_pair
    cdef Junction * jj
    cdef pair[string, Exon *] ex_pair
    cdef Exon * ex
    cdef string gne_id = gne.get_id()

    sg_gene(db, gne_id, gne.get_name(), string(1, gne.get_strand()), gne.get_chromosome())
    for jj_pair in gne.junc_map_:

        jj = jj_pair.second
        if not jj.get_denovo_bl(): continue
        if jj.get_start() == C_FIRST_LAST_JUNC:
            sg_alt_start(db, gne_id, jj.get_end())
            continue

        if jj.get_end() == C_FIRST_LAST_JUNC:
            sg_alt_end(db, gne_id, jj.get_start())
            continue
        # with gil:
        #     print("## ", gne_id, jj.get_start(), jj.get_end(), jj.get_annot(), jj.get_simpl_fltr())
        sg_junction(db, gne_id, jj.get_start(), jj.get_end(), jj.get_annot(), jj.get_simpl_fltr(), jj.get_constitutive())

    for ex_pair in gne.exon_map_:
        ex = ex_pair.second
        # with gil:
        #     print(ex_pair.first, ex.get_start(), ex.get_end())
        sg_exon(db, gne_id, ex.get_start(), ex.get_end(), ex.db_start_, ex.db_end_, ex.annot_ )
        if ex.has_out_intron():
            ir = ex.ob_irptr
            if ir.get_ir_flag():
                # with gil:
                #     print(gne_id, ir.get_start(), ir.get_end(), ir.get_annot(), ir.is_connected())
                sg_intron_retention(db, gne_id, ir.get_start(), ir.get_end(), ir.get_annot(), ir.get_simpl_fltr(),
                                    ir.get_constitutive())


cdef int simplify(list file_list, map[string, Gene*] gene_map, vector[string] gid_vec,
                   map[string, overGene_vect_t] gene_list, object conf, object logger) except -1 :

    cdef int nsamples = len(file_list)
    cdef int nthreads = conf.nthreads
    cdef int denovo_simpl = conf.simpl_denovo
    cdef int db_simple = conf.simpl_db
    cdef int ir_simpl = conf.simpl_ir
    cdef int i, j
    cdef int n = gene_map.size()
    cdef int njunc
    cdef int strandness
    cdef unsigned int irbool, coord1, coord2, sreads, npos
    cdef unsigned int min_experiments
    cdef map[string, int] junc_tlb
    cdef Gene * gg
    cdef Gene_vect_t gene_l
    cdef np.float32_t simpl_fltr = conf.simpl_psi
    cdef bint irb =  conf.ir
    cdef bint val = True ;
    cdef bint bsimpl = (conf.simpl_psi >= 0)
    cdef int last_it_grp;
    cdef Intron * ir_ptr;
    cdef string chrom, gid, key, jid
    cdef char strand

    logger.info('Starting simplification %s' % bsimpl)
    for tmp_str, group_list in conf.tissue_repl.items():
        name = tmp_str.encode('utf-8')
        last_it_grp = group_list[len(group_list) - 1]
        min_experiments = conf.min_experiments[tmp_str]

        for i in group_list:

            strandness = conf.strand_specific[file_list[i][0]]
            if file_list[i][2]:
                junc_file =  file_list[i][1]
            else:
                junc_file = "%s/%s.%s" % (conf.outDir, file_list[i][0], JUNC_FILE_FORMAT)

            with open(junc_file, 'rb') as fp:
                junc_ids = np.load(fp)['junc_info']
                njunc = junc_ids.shape[0]

                for j in prange(njunc, nogil=True, num_threads=nthreads):
                    gene_l = Gene_vect_t()
                    with gil:
                        jid     = junc_ids[j][0]
                        coord1  = junc_ids[j][1]
                        coord2  = junc_ids[j][2]
                        sreads  = junc_ids[j][3]
                        irbool  = junc_ids[j][5]
                        chrom   = jid.split(b':')[0]
                        strand  = <char> jid.split(b':')[1][0]

                    if irbool == 0:
                        with gil:
                            junc_tlb[jid] = sreads

                    elif irb:
                        find_gene_from_junc(gene_list, chrom, strand, coord1, coord2, gene_l, irbool, bsimpl)
                        with gil:
                            gid = b':'.join(jid.split(b':')[3:])
                        for gg in gene_l:
                            if gg.get_id() != gid:
                                continue
                            irv = find_intron_retention(gg, coord1, coord2)
                            for ir_ptr in irv:
                                key = key_format(gg.get_id(), ir_ptr.get_start(), ir_ptr.get_end(), val)
                                with gil:
                                    junc_tlb[key] = sreads
                    gene_l.clear()
                del junc_ids
            logger.debug('Simplifying file %s %s %s %s' %(file_list[i][0], (i==last_it_grp), last_it_grp, i))
            for j in prange(n, nogil=True, num_threads=nthreads):
                gg = gene_map[gid_vec[j]]
                gg.simplify(junc_tlb, simpl_fltr, strandness, denovo_simpl, db_simple, ir_simpl,
                            (i==last_it_grp), min_experiments)
            junc_tlb.clear()
    logger.info('Finished simplification')



## OPEN API FOR PYTHON
cdef _core_build(str transcripts, list file_list, object conf, object logger):
    cdef int n, i
    cdef int nthreads = conf.nthreads
    # cdef map[string, vector[string]] gene_junc_tlb
    cdef map[string, int] lsv_juncs_tlb
    cdef vector[LSV*] out_lsvlist
    cdef int nsamples = len(file_list)
    cdef int k=conf.k, m=conf.m
    cdef bint ir = conf.ir
    cdef int nlsv
    cdef map[string, Gene*] gene_map
    cdef map[string, overGene_vect_t] gene_list

    cdef vector[string] gid_vec
    cdef string sg_filename = get_builder_splicegraph_filename(conf.outDir).encode('utf-8')
    cdef string fname
    cdef string outDir = conf.outDir.encode('utf-8')
    cdef int strandness, cnt
    cdef sqlite3* db
    cdef bint bsimpl = (conf.simpl_psi >= 0)
    cdef bint dumpCJunctions = conf.dump_const_j
    cdef vector[string] cjuncs
    cdef bint enable_anot_ir = conf.annot_ir_always

    logger.info("Parsing GFF3")
    majiq_io.read_gff(transcripts, gene_map, gid_vec, bsimpl, enable_anot_ir, logger)

    prepare_genelist(gene_map, gene_list)
    n = gene_map.size()
    init_splicegraph(sg_filename, conf)
    logger.info("Reading bamfiles")
    _find_junctions(file_list, gene_map, gid_vec, gene_list, conf, logger)

    if conf.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss)/(1024**2)
        logger.info("PRE LOOP Memory used %.2f MB" % mem_allocated)

    if conf.juncfiles_only:
        return
    logger.info("Detecting LSVs ngenes: %s " % n)
    open_db(sg_filename, &db)

    for i in prange(n, nogil=True, num_threads=nthreads):
        gg = gene_map[gid_vec[i]]
        with gil:
            logger.debug("%s] Detect exons" % gg.get_id())
        gg.detect_exons()
        if ir:
            with gil:
                logger.debug("%s] Connect introns" % gg.get_id())
            gg.connect_introns()

    if bsimpl:
        simplify(file_list, gene_map, gid_vec, gene_list, conf, logger)
    for i in prange(n, nogil=True, num_threads=nthreads):
        gg = gene_map[gid_vec[i]]
        gg.get_constitutive_junctions(cjuncs)
        gene_to_splicegraph(gg, db)
        with gil:
            logger.debug("[%s] Detect LSVs" % gg.get_id())
        nlsv = gg.detect_lsvs(out_lsvlist)

    if cjuncs.size()>0 and dumpCJunctions:
        with open("%s/constitutive_junctions.tsv" % conf.outDir, 'w+') as fp:
            fp.write('#GENEID\tCHROMOSOME\tJUNC_START\tJUNC_END\tDONOR_START\tDONOR_END\tACCEPTOR_START\tACCEPTOR_END\n')
            for jstr in cjuncs:
                fp.write("%s\n" % jstr.decode('utf8'))


    logger.debug("Generate TLB")
    fill_junc_tlb(out_lsvlist, lsv_juncs_tlb)

    logger.info("%s LSV found" % out_lsvlist.size())
    if conf.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss)/(1024**2)
        logger.info("POST LOOP Memory used %.2f MB" % mem_allocated)

    for i in range(nsamples):
        strandness = conf.strand_specific[file_list[i][0]]
        if conf.mem_profile:
            mem_allocated = int(psutil.Process().memory_info().rss)/(1024**2)
            logger.info("PRE OUT Memory used %.2f MB" % mem_allocated)
        cnt  = _output_majiq_file(out_lsvlist, gene_list, lsv_juncs_tlb, file_list[i], outDir, db, m, ir, bsimpl,
                                  logger, nthreads)

        logger.info('%s: %d LSVs' %(file_list[i][0], cnt))
        if conf.mem_profile:
            mem_allocated = int(psutil.Process().memory_info().rss)/(1024**2)
            logger.info("POST OUT  Memory used %.2f MB" % mem_allocated)

    close_db(db)
    free_genelist(gene_list)
    free_lsvlist(out_lsvlist)


def build(args):
    pipeline_run(Builder(args))

class Builder(BasicPipeline):

    def run(self):
        # if self.simplify is not None and len(self.simplify) not in (0, 2):
        #     raise RuntimeError('Simplify requires 2 values type of junctions afected and E(PSI) threshold.')
        if not os.path.exists(self.conf):
            raise RuntimeError("Config file %s does not exist" % self.conf)
        majiq_config = Config(self.conf, self)
        self.builder(majiq_config)

    def builder(self, majiq_config):

        logger = majiq_logger.get_logger("%s/majiq.log" % majiq_config.outDir, silent=self.silent, debug=self.debug)
        logger.info("Majiq Build v%s-%s" % (VERSION, get_git_version()))
        logger.info("Command: %s" % " ".join(sys.argv))

        _core_build(self.transcripts, majiq_config.sam_list, majiq_config, logger)

        if self.mem_profile:
            mem_allocated = int(psutil.Process().memory_info().rss)/(1024**2)
            logger.info("Max Memory used %.2f MB" % mem_allocated)
        logger.info("MAJIQ Builder is ended successfully!")
