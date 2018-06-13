from majiq.src.internals.io_bam cimport IOBam, gene_vect_t
from majiq.src.internals.grimoire cimport Junction, Gene, Exon, LSV, Jinfo, Intron
from majiq.src.internals.grimoire cimport sortGeneList, find_intron_retention
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

from libc.stdlib cimport calloc, free
from numpy.random  import choice
from scipy.stats import nbinom, poisson
import numpy as np
cimport numpy as np

cdef _store_junc_file(np.ndarray boots, list junc_ids, str experiment_name, str outDir):

    cdef str out_file = "%s/%s.juncs" % (outDir, experiment_name)
    cdef dict vals = {'bootstrap': boots}
    dt = np.dtype('S250, u4, u4, f4, f4, u4')
    vals['junc_info'] = np.array(junc_ids, dtype=dt)
    with open(out_file, 'w+b') as ofp:
        np.savez(ofp, **vals)

def update_splicegraph_junction(sg, gene_id, start, end, nreads, exp):
    if FIRST_LAST_JUNC not in (start, end):
        sg.junction(gene_id, start, end).update_reads(exp, nreads)

cdef int _output_lsv_file_single(vector[LSV*] out_lsvlist, str experiment_name, map[string, vector[string]] tlb_j_g,
                                 map[string, gene_vect_t] gene_map, str outDir, int nthreads, unsigned int msamples,
                                 bint irb, int strandness):
    cdef dict cov_dict = {}
    cdef int nlsv = out_lsvlist.size()
    cdef str out_file, junc_file
    cdef list type_list = []
    cdef LSV* lsv_ptr
    cdef int njunc
    cdef Jinfo jobj_ptr
    cdef map[string, Jinfo] tlb_juncs
    cdef map[string, Jinfo] tlb_ir
    cdef int i, j, junc_idx
    cdef list junc_info = []
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] boots
    cdef np.ndarray junc_ids
    cdef object all_juncs
    cdef string key
    cdef string lsvid
    cdef string ir_key
    cdef vector[Intron *] irv
    cdef Intron * ir_ptr

    junc_file = "%s/%s.juncs" % (outDir, experiment_name)
    #print("###", junc_file, get_builder_splicegraph_filename(outDir))
    out_file = "%s/%s.majiq" % (outDir, experiment_name)

    with SpliceGraph(get_builder_splicegraph_filename(outDir)) as sg:
        with open(junc_file, 'rb') as fp:
            all_juncs = np.load(fp)
            boots = all_juncs['bootstrap']
            junc_ids = all_juncs['junc_info']

            ''' If we move towards this solution we can remove elements in Jinfo and make it smaller'''
            #TODO: with tlb_j_g we can filter junctions are no present in this case like smaller annot db
            for i in range(junc_ids.shape[0]):

                if junc_ids[i][5] == 0:
                    # print("###: %s" , junc_ids[i][0])
                    tlb_juncs[junc_ids[i][0]] = Jinfo(i, junc_ids[i][1], junc_ids[i][2], junc_ids[i][3], junc_ids[i][4])
                    if tlb_j_g.count(junc_ids[i][0]) == 0 :
                        continue

                    for gid in tlb_j_g[junc_ids[i][0]]:
                        # print("ADD COUNT", i, gid.decode('utf-8'), junc_ids[i][1], junc_ids[i][2], junc_ids[i][3], junc_ids[i][0], experiment_name)
                        update_splicegraph_junction(sg, gid.decode('utf-8'), junc_ids[i][1], junc_ids[i][2],
                                                    junc_ids[i][3], experiment_name)


                elif irb:
                    #intron retention case
                    # print('KKK2')
                    key = junc_ids[i][0].split(b':')[0]
                    irv = find_intron_retention(gene_map[key], junc_ids[i][1], junc_ids[i][2])

                    kk = [ print(ir_ptr.get_gene(), ir_ptr.get_start(), ir_ptr.get_end()) for ir_ptr in irv ]

                    print("######\n")
                    for ir_ptr in irv:

                        sg.intron_retention(ir_ptr.get_gene().get_id().decode('utf-8'), ir_ptr.get_start(),
                                            ir_ptr.get_end()).update_reads(experiment_name, junc_ids[i][3])
                        tlb_ir[ir_ptr.get_key(ir_ptr.get_gene())] = Jinfo(i, junc_ids[i][1], junc_ids[i][2],
                                                                  junc_ids[i][3], junc_ids[i][4])

            del junc_ids

            for j in prange(nlsv, nogil=True, num_threads=nthreads):
                lsv_ptr = out_lsvlist[j]
                njunc = lsv_ptr.get_num_variations()
                lsvid = lsv_ptr.get_id()
                with gil:
                    cov_dict[lsvid.decode('utf-8')] = np.zeros(shape=(njunc, msamples), dtype=np.float32)
                    type_list.append((lsvid.decode('utf-8'), lsv_ptr.get_type()))
                junc_idx = 0

                for junc in lsv_ptr.get_junctions():
                    key = junc.get_key(lsv_ptr.get_gene(), strandness)
                    # with gil:
                        # print(key)
                    if tlb_juncs.count(key) > 0 :
                        jobj_ptr = tlb_juncs[key]
                        with gil:
                            # print (boots[tlb_juncs[key].index])
                            cov_dict[lsvid.decode('utf-8')][junc_idx] = boots[tlb_juncs[key].index]
                            junc_info.append((lsvid.decode('utf-8'), junc.get_start(), junc.get_end(),
                                              jobj_ptr.sreads, jobj_ptr.npos))
                    junc_idx = junc_idx + 1

                ir_ptr = lsv_ptr.get_intron()
                if irb and ir_ptr != <Intron * > 0:
                    key = ir_ptr.get_key(lsv_ptr.get_gene())
                    if tlb_ir.count(key) > 0 :
                        jobj_ptr = tlb_ir[key]
                        with gil:
                            cov_dict[lsvid.decode('utf-8')][junc_idx] = boots[jobj_ptr.index]
                            junc_info.append((lsvid.decode('utf-8'), ir_ptr.get_start(), ir_ptr.get_end(),
                                              jobj_ptr.sreads, jobj_ptr.npos))

            majiq_io.dump_lsv_coverage(out_file, cov_dict, type_list, junc_info, experiment_name)
            tlb_juncs.clear()
            tlb_ir.clear()

    return len(type_list)


cdef _find_junctions(list file_list, vector[Gene*] gene_vec,  object conf, object logger):

    cdef int n = gene_vec.size()
    cdef int nthreads = conf.nthreads
    cdef int nsamples = len(file_list)
    cdef unsigned int minpos = conf.minpos
    cdef unsigned int minreads = conf.minreads
    cdef unsigned int denovo_thresh= conf.min_denovo
    cdef float min_ir_cov = conf.min_intronic_cov
    cdef bint ir = conf.ir
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
    cdef unsigned int jlimit

    cdef int* jvec

    cdef pair[string, gene_vect_t] vector_gene
    cdef map[string, gene_vect_t] gene_list
    cdef map[string, unsigned int] j_ids
    cdef map[string, int*] junc_summary
    cdef pair[string, unsigned int] it
    cdef map[string, vector[string]] gene_junc_tlb

    # print("SIZE", gene_vec.size())
    sortGeneList(gene_vec)

    for gobj in gene_vec:
        if gene_list.count(gobj.get_chromosome()) == 0:
            gene_list[gobj.get_chromosome()] = gene_vect_t()
        gene_list[gobj.get_chromosome()].push_back(gobj)

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
                    c_iobam.detect_introns(min_ir_cov, min_experiments, 0.8)
                njunc = c_iobam.get_njuncs()
                with gil:
                        logger.debug('Total Junctions and introns %s' %(njunc))

            boots = np.zeros(shape=(njunc, m), dtype=np.float32)
            with nogil:
                c_iobam.boostrap_samples(m, k, <np.float32_t *> boots.data)
                j_ids = c_iobam.get_junc_map()
                jvec = c_iobam.get_junc_vec_summary()
                jlimit = c_iobam.get_junc_limit_index()

            logger.debug("Update flags")
            for i in prange(n, nogil=True, num_threads=nthreads):
                gg = gene_vec[i]
                gg.update_junc_flags(eff_len, (j==last_it_grp), minreads, minpos, denovo_thresh, min_experiments)

            logger.debug("Done Update flags")
            junc_ids = [0] * njunc

            jv_idx = 0
            for it in j_ids:

                # print(file_list[j][0], j, last_it_grp, (j==last_it_grp), it.second, it.first)
                tmp_str = it.first.decode('utf-8').split(':')[2]
                start, end = (int(xx) for xx in tmp_str.split('-'))
                junc_ids[it.second] = (it.first.decode('utf-8'), start, end, jvec[it.second], jvec[it.second + njunc], int(it.second>= jlimit))

            logger.info('Done Reading file %s' %(file_list[j][0]))
            _store_junc_file(boots, junc_ids, file_list[j][0], conf.outDir)


    logger.info("Detecting LSVs ngenes: %s " % n)
    for i in prange(n, nogil=True, num_threads=nthreads):
        gg = gene_vec[i]
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
        # with gil:
        #     gene_to_splicegraph(gg, conf)
        with gil:
            logger.debug("[%s] Detect LSVs" % gg.get_id())
        nlsv = gg.detect_lsvs(out_lsvlist)


    logger.info("Generating Splicegraph")
    for i in range(n):
        gene_to_splicegraph(gene_vec[i], conf)

    logger.info("%s LSV found" % out_lsvlist.size())
    for j in range(nsamples):
        strandness = conf.strand_specific[file_list[j][0]]
        cnt = _output_lsv_file_single(out_lsvlist, file_list[j][0], gene_junc_tlb, gene_list, conf.outDir, nthreads, m, ir, strandness)
        logger.info('%s: %d LSVs' %(file_list[j][0], cnt))



from voila.api import SpliceGraph
from majiq.src.config import Config

cdef init_splicegraph(filename):
    majiq_config = Config()

    # erase splice graph file
    with SpliceGraph(filename, delete=True) as sg:
        print([(xx, type(xx)) for xx in majiq_config.exp_list])
        sg.experiment_names = majiq_config.exp_list
        sg.genome = majiq_config.genome

cdef gene_to_splicegraph(Gene * gne, majiq_config):

    cdef pair[string, Junction *] jj_pair
    cdef Junction * jj
    cdef pair[string, Exon *] ex_pair
    cdef Exon * ex
    cdef str gne_id

    alt_empty_starts = []
    alt_empty_ends = []
    gne_id = gne.get_id().decode('utf-8')
    # print('GENEID:', gne_id, gne.get_name(), chr(gne.get_strand()), gne.get_chromosome())
    with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir)) as sg:
        sg.gene(gne_id).add(name=gne.get_name().decode('utf-8'), strand=chr(gne.get_strand()), chromosome=gne.get_chromosome().decode('utf-8'))

        for jj_pair in gne.junc_map_:
            jj = jj_pair.second
            # if gne_id == 'ENSMUSG00000006498':
            #     print(jj_pair.first, jj.get_denovo_bl(), jj.get_start(), jj.get_end(), jj.get_annot())

            if not jj.get_denovo_bl(): continue
            if jj.get_start() == FIRST_LAST_JUNC:
                alt_empty_starts.append(jj.get_end())
                continue

            if jj.get_end() == FIRST_LAST_JUNC:
                alt_empty_ends.append(jj.get_start())
                continue

            # TODO: add transcripts
            # if gne_id == 'ENSMUSG00000006498':
            #     print('JUNC', gne_id, jj.get_start(), jj.get_end())
            sg.junction(gne_id, jj.get_start(), jj.get_end()).add(annotated=jj.get_annot())

        for ex_pair in gne.exon_map_:
            ex = ex_pair.second
        #for ex in sorted(exon_dict, key=lambda x: (x.start, x.end)):
            # if ex.intron_:
            #     sg.exon(gne_id, ex.start_, ex.end_).add(annotated=ex.annot_, intron_retention=True)
            #     continue

            alt_start = []
            for jj in ex.ib:
                if jj.get_end() in alt_empty_starts:
                    alt_start.append(jj.get_end())

            alt_ends = []
            for jj in ex.ob:
                if jj.get_start() in alt_empty_ends:
                    alt_ends.append(jj.get_start())

            extra_coords = []
            if ex.annot_:
                if ex.get_start() < ex.db_start_:
                    extra_coords.append([ex.get_start(), ex.db_start_ - 1])
                if ex.get_end() > ex.db_end_:
                    extra_coords.append([ex.db_end_ + 1, ex.get_end()])

            sg.exon(gne_id, ex.get_start(), ex.get_end()).add(coords_extra=extra_coords, annotated=ex.annot_,
                                                    alt_starts=alt_start, alt_ends=alt_ends)

        for ir in gne.intron_vec_:
            # print("INTRON", gne_id,  ir.get_start(), ir.get_end())
            sg.intron_retention(gne_id, ir.get_start(), ir.get_end()).add(annotated=ir.get_annot())




## OPEN API FOR PYTHON
cdef _core_build(str transcripts, list file_list, object conf, object logger):
    cdef vector[Gene *] gene_vec
    majiq_io.read_gff(transcripts, gene_vec, logger)
    init_splicegraph(get_builder_splicegraph_filename(conf.outDir))
    _find_junctions(file_list, gene_vec, conf, logger)


cpdef core_build(str transcripts, list file_list, object conf, object logger):
    _core_build(transcripts, file_list, conf, logger)





