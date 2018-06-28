from majiq.src.internals.io_bam cimport IOBam, gene_vect_t
from majiq.src.internals.grimoire cimport Junction, Gene, Exon, LSV, Jinfo, Intron
from majiq.src.internals.grimoire cimport sortGeneList, find_intron_retention
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

# from majiq.src.polyfitnb import fit_nb

from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from cython.parallel import prange



cdef extern from "sqlite3.h":
    struct sqlite3

import numpy as np
cimport numpy as np

cdef int C_FIRST_LAST_JUNC = FIRST_LAST_JUNC

cdef _store_junc_file(np.ndarray boots, list junc_ids, str experiment_name, str outDir):

    cdef str out_file = "%s/%s.juncs" % (outDir, experiment_name)
    cdef dict vals = {'bootstrap': boots}
    dt = np.dtype('S250, u4, u4, f4, f4, u4')
    vals['junc_info'] = np.array(junc_ids, dtype=dt)
    with open(out_file, 'w+b') as ofp:
        np.savez(ofp, **vals)

cdef void update_splicegraph_junction(sqlite3 *db, string gene_id, int start, int end, int nreads, string exp) nogil:
    if C_FIRST_LAST_JUNC != start and C_FIRST_LAST_JUNC != end:
        sg_junction_reads(db, nreads, exp, gene_id, start, end)

cdef int _output_lsv_file_single(vector[LSV*] out_lsvlist, string experiment_name, map[string, vector[string]] tlb_j_g,
                                 map[string, Gene*] gene_map, string outDir, int nthreads, unsigned int msamples,
                                 bint irb, int strandness, object logger) :

    # cdef np.float32_t[:,:,:] boots
    # cdef np.float32_t* boots
    cdef unsigned int irbool, coord1, coord2, sreads, npos
    cdef string jid
    cdef string geneid
    cdef Gene* gneObj
    cdef vector[Intron *] irv
    cdef Intron * ir_ptr
    cdef sqlite3* db
    cdef string sg_filename
    cdef LSV* lsv_ptr
    cdef string lsvid, gid
    cdef unsigned int i, j, junc_idx

    cdef int nlsv = out_lsvlist.size()

    cdef int njunc
    cdef Jinfo jobj_ptr
    cdef map[string, Jinfo] tlb_juncs
    cdef map[string, Jinfo] tlb_ir


    cdef string key
    cdef string ir_key
    cdef char strand
    cdef Junction * jObj

    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
    cdef np.ndarray junc_ids
    cdef str out_file, junc_file
    cdef dict cov_dict = {}
    cdef list type_list = []
    cdef list junc_info = []
    cdef object all_juncs

    # cdef map[string, np.float32_t*] cov_dict

    # with gil:

    sg_filename = get_builder_splicegraph_filename(outDir).encode('utf-8')
    junc_file = "%s/%s.juncs" % (outDir.decode('utf-8'), experiment_name.decode('utf-8'))
    out_file = "%s/%s.majiq" % (outDir.decode('utf-8'), experiment_name.decode('utf-8'))
    with open(junc_file, 'rb') as fp:
        all_juncs = np.load(fp)
        boots = all_juncs['bootstrap']
        junc_ids = all_juncs['junc_info']
        njunc = junc_ids.shape[0]

    cov_dict = {}
    junc_info = []
    type_list = []

    with nogil:
        db = open_db(sg_filename)
        for i in range(njunc):
            with gil:
                jid     = junc_ids[i][0]
                coord1  = junc_ids[i][1]
                coord2  = junc_ids[i][2]
                sreads  = junc_ids[i][3]
                npos    = junc_ids[i][4]
                irbool  = junc_ids[i][5]

            jobj_ptr = Jinfo(i, coord1, coord2, sreads, npos)

            if irbool == 0:
                tlb_juncs[jid] = jobj_ptr
                if tlb_j_g.count(jid) > 0:
                    for gid in tlb_j_g[jid]:
                        # print("ADD COUNT", i, gid.decode('utf-8'), junc_ids[i][1], junc_ids[i][2], junc_ids[i][3],
                        # junc_ids[i][0], experiment_name)
                        update_splicegraph_junction(db, gid, coord1, coord2, sreads, experiment_name)

            elif irb:
                # print("###: %s" , junc_ids[i][0])
                with gil:
                    geneid = jid.split(b':')[3]
                if gene_map.count(geneid) > 0 :
                    gneObj = gene_map[geneid]
                    gid = gneObj.get_id()
                    irv = find_intron_retention(gneObj, coord1, coord2)

                    for ir_ptr in irv:
                        sg_intron_retention_reads(db, sreads, experiment_name, geneid, ir_ptr.get_start(),
                                                  ir_ptr.get_end())
                        tlb_ir[ir_ptr.get_key(ir_ptr.get_gene())] = jobj_ptr

        close_db(db)
    # with gil:
    del junc_ids

    logger.info("Create majiq file")

    # for j in prange(nlsv, nogil=True, num_threads=nthreads):
    with nogil:
        for j in range(nlsv):
            lsv_ptr = out_lsvlist[j]
            njunc = lsv_ptr.get_num_variations()
            lsvid = lsv_ptr.get_id()
            with gil:
                # cov_dict[lsvid] = np.zeros(shape=(njunc, msamples), dtype=np.float32)
                type_list.append((lsvid.decode('utf-8'), lsv_ptr.get_type()))
            junc_idx = 0

            for junc in lsv_ptr.get_junctions():
                key = junc.get_key(lsv_ptr.get_gene(), strandness)
                # with gil:
                    # print(key)
                if tlb_juncs.count(key) > 0 :
                    jobj_ptr = tlb_juncs[key]
                    sreads = jobj_ptr.sreads
                    npos = jobj_ptr.npos
                else:
                    sreads = 0
                    npos = 0
                with gil:
                        # print (boots[tlb_juncs[key].index])
                    cov_dict[lsvid][junc_idx] = boots[tlb_juncs[key].index]
                    junc_info.append((lsvid.decode('utf-8'), junc.get_start(), junc.get_end(),
                                          jobj_ptr.sreads, jobj_ptr.npos))
                junc_idx = junc_idx + 1

            ir_ptr = lsv_ptr.get_intron()
            if irb and ir_ptr != <Intron * > 0:
                key = ir_ptr.get_key(lsv_ptr.get_gene())
                if tlb_ir.count(key) > 0 :
                    jobj_ptr = tlb_ir[key]
                    sreads = jobj_ptr.sreads
                    npos = jobj_ptr.npos
                else:
                    sreads = 0
                    npos = 0

                with gil:
                    cov_dict[lsvid][junc_idx] = boots[jobj_ptr.index]
                    junc_info.append((lsvid.decode('utf-8'), ir_ptr.get_start(), ir_ptr.get_end(),
                                      jobj_ptr.sreads, jobj_ptr.npos))
    # with gil:
    logger.info("Dump majiq file")
    majiq_io.dump_lsv_coverage(out_file, cov_dict, type_list, junc_info, experiment_name.decode('utf-8'))
    nlsv = len(type_list)

    tlb_juncs.clear()
    tlb_ir.clear()

    return nlsv


cdef _find_junctions(list file_list, map[string, Gene*]& gene_map, vector[string] gid_vec, object conf, object logger):

    cdef int n = gene_map.size()
    cdef int nthreads = conf.nthreads

    cdef unsigned int minpos = conf.minpos
    cdef unsigned int minreads = conf.minreads
    cdef unsigned int denovo_thresh= conf.min_denovo
    cdef float min_ir_cov = conf.min_intronic_cov
    cdef int k=conf.k, m=conf.m
    cdef float pvalue_limit=conf.pvalue_limit
    cdef unsigned int min_experiments = 1 if conf.min_exp == -1 else conf.min_exp
    cdef unsigned int eff_len = conf.readLen - 2*MIN_BP_OVERLAP
    cdef bint ir = conf.ir

    cdef int i, j
    cdef int strandness, njunc
    cdef list group_list
    cdef int last_it_grp

    cdef Gene * gg
    cdef IOBam c_iobam
    cdef string name, bamfile
    cdef str tmp_str
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
    cdef list junc_ids
    cdef float fitfunc_r
    cdef unsigned int jlimit

    cdef int* jvec

    cdef pair[string, gene_vect_t] vector_gene
    cdef map[string, gene_vect_t] gene_list
    cdef map[string, unsigned int] j_ids
    cdef pair[string, unsigned int] it
    cdef pair[string, Gene *] git
    # sortGeneList(gene_vec)

    for git in gene_map:
        gg = git.second
        if gene_list.count(gg.get_chromosome()) == 0:
            gene_list[gg.get_chromosome()] = gene_vect_t()
        gene_list[gg.get_chromosome()].push_back(gg)

    for vector_gene in gene_list:
        sortGeneList(vector_gene.second)

    for tmp_str, group_list in conf.tissue_repl.items():
        name = tmp_str.encode('utf-8')
        last_it_grp = group_list[len(group_list) - 1]
        for j in group_list:
            logger.info('Reading file %s' %(file_list[j][0]))
            bamfile = ('%s/%s.%s' % (conf.sam_dir, file_list[j][0], SEQ_FILE_FORMAT)).encode('utf-8')
            strandness = conf.strand_specific[file_list[j][0]]

            with nogil:
                c_iobam = IOBam(bamfile, strandness, eff_len, nthreads, gene_list)
                c_iobam.ParseJunctionsFromFile(False)
                if ir:
                    with gil:
                        logger.info('Detect Intron retention %s' %(file_list[j][0]))
                    c_iobam.detect_introns(min_ir_cov, min_experiments, 0.8, (j==last_it_grp))
                njunc = c_iobam.get_njuncs()
                with gil:
                        logger.debug('Total Junctions and introns %s' %(njunc))

            # fitfunc_r = fit_nb(c_iobam.junc_vec, logger=logger)
            fitfunc_r = 0
            boots = np.zeros(shape=(njunc, m), dtype=np.float32)
            with nogil:
                c_iobam.boostrap_samples(m, k, <np.float32_t *> boots.data, fitfunc_r, pvalue_limit)
                j_ids = c_iobam.get_junc_map()
                jvec = c_iobam.get_junc_vec_summary()
                jlimit = c_iobam.get_junc_limit_index()

            logger.debug("Update flags")
            for i in prange(n, nogil=True, num_threads=nthreads):
                gg = gene_map[gid_vec[i]]
                gg.update_junc_flags(eff_len, (j==last_it_grp), minreads, minpos, denovo_thresh, min_experiments)
            logger.debug("Done Update flags")
            junc_ids = [0] * njunc
            for it in j_ids:
                # print(file_list[j][0], j, last_it_grp, (j==last_it_grp), it.second, it.first)
                tmp_str = it.first.decode('utf-8').split(':')[2]
                start, end = (int(xx) for xx in tmp_str.split('-'))
                junc_ids[it.second] = (it.first.decode('utf-8'), start, end, jvec[it.second], jvec[it.second + njunc],
                                       int(it.second>= jlimit))

            logger.info('Done Reading file %s' %(file_list[j][0]))
            _store_junc_file(boots, junc_ids, file_list[j][0], conf.outDir)
            c_iobam.free_iobam()




cdef init_splicegraph(string filename, object conf):

    # erase splice graph file
    print(filename)
    with SpliceGraph(filename.decode('utf-8'), delete=True) as sg:
        sg.experiment_names = conf.exp_list
        sg.genome = conf.genome

cdef void gene_to_splicegraph(Gene * gne, string sg_filename) nogil:

    cdef pair[string, Junction *] jj_pair
    cdef Junction * jj
    cdef pair[string, Exon *] ex_pair
    cdef Exon * ex
    cdef string gne_id
    cdef sqlite3* db

    gne_id = gne.get_id()

    db = open_db(sg_filename)
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
        sg_junction(db, gne_id, jj.get_start(), jj.get_end(), jj.get_annot())

    for ex_pair in gne.exon_map_:
        ex = ex_pair.second
        sg_exon(db, gne_id, ex.get_start(), ex.get_end(), ex.db_start_, ex.db_end_, ex.annot_ )

    for ir in gne.intron_vec_:
        if not ir.get_ir_flag(): continue
        sg_intron_retention(db, gne_id, ir.get_start(), ir.get_end(), ir.get_annot())
    close_db(db)


## OPEN API FOR PYTHON
cdef _core_build(str transcripts, list file_list, object conf, object logger):
    # cdef vector[Gene *] gene_vec
    cdef int n, i
    cdef int nthreads = conf.nthreads
    cdef map[string, vector[string]] gene_junc_tlb
    cdef vector[LSV*] out_lsvlist
    cdef int nsamples = len(file_list)
    cdef int k=conf.k, m=conf.m
    cdef bint ir = conf.ir
    cdef int nlsv
    cdef map[string, Gene*] gene_map
    cdef vector[string] gid_vec
    cdef string sg_filename = get_builder_splicegraph_filename(conf.outDir).encode('utf-8')
    cdef string fname
    cdef string outDir = conf.outDir.encode('utf-8')
    cdef int strandness, cnt

    majiq_io.read_gff(transcripts, gene_map, gid_vec, logger)
    # logger.info("PRE INIT SG1")
    n = gene_map.size()
    # logger.info("PRE INIT SG2 %s" % n)
    init_splicegraph(sg_filename, conf)
    # logger.info("POST INIT SG")
    _find_junctions(file_list, gene_map, gid_vec, conf, logger)

    logger.info("Detecting LSVs ngenes: %s " % n)
    for i in prange(n, nogil=True, num_threads=nthreads):
        gg = gene_map[gid_vec[i]]
        with gil:
            logger.debug("%s] Detect exons" % gg.get_id())
        gg.detect_exons()
        if ir:
            with gil:
                logger.debug("%s] Connect introns" % gg.get_id())
            gg.connect_introns()

        with gil:
            logger.debug("[%s] Generate TLB" % gg.get_id())

        gg.fill_junc_tlb(gene_junc_tlb)
        gene_to_splicegraph(gg, sg_filename)
        # with gil:
        #     gene_to_splicegraph2(gene_map[gid_vec[i]], sg_filename)
        with gil:
            logger.debug("[%s] Detect LSVs" % gg.get_id())

        nlsv = gg.detect_lsvs(out_lsvlist)

    logger.info("%s LSV found" % out_lsvlist.size())
    for i in prange(nsamples, nogil=True, num_threads=nthreads):
        with gil:
            fname = file_list[i][0].encode('utf-8')
            strandness = conf.strand_specific[file_list[i][0]]
            cnt = _output_lsv_file_single(out_lsvlist, fname, gene_junc_tlb, gene_map, outDir,
                                          nthreads, m, ir, strandness, logger)
            logger.info('%s: %d LSVs' %(fname.decode('utf-8'), cnt))

cpdef core_build(str transcripts, list file_list, object conf, object logger):
    _core_build(transcripts, file_list, conf, logger)





