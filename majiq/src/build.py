import multiprocessing as mp
import os
import sys

import h5py

import majiq.src.io as majiq_io
from majiq.src.io_bam import read_juncs, find_introns
import majiq.src.io_bam as majiq_io_bam

from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper

import majiq.src.logger as majiq_logger
import majiq.src.multiproc as majiq_multi
import majiq.src.normalize as majiq_norm
from majiq.grimoire.exon import detect_exons
from majiq.grimoire.junction import create_junction
from majiq.grimoire.lsv import detect_lsvs
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.config import Config
from majiq.src.constants import *
from majiq.src.polyfitnb import fit_nb
from majiq.src.voila_wrapper import gene_to_splicegraph, init_splicegraph, update_splicegraph_junctions
import math
import datetime


def build(args):
    pipeline_run(Builder(args))


def find_new_junctions(file_list, chunk, process_conf, logger):

    majiq_config = Config()
    dict_of_genes = majiq_io.retrieve_db_genes(majiq_config.outDir)
    list_exons = {}
    dict_junctions = {}
    for gne_id, gene_obj in dict_of_genes.items():
        list_exons[gne_id] = []
        dict_junctions[gne_id] = {}
        majiq_io.retrieve_db_info(gne_id, majiq_config.outDir, list_exons[gne_id],
                                  dict_junctions[gne_id], None)
        detect_exons(dict_junctions[gne_id], list_exons[gne_id])

    for is_junc_file, fname, name in file_list:
        logger.info('READ JUNCS from %s, %s' % (fname, majiq_config.strand_specific))
        read_juncs(fname, is_junc_file, list_exons, dict_of_genes, dict_junctions,
                   majiq_config.strand_specific, process_conf.queue, gname=name)


def find_new_introns(file_list, chunk, process_conf, logger):

    majiq_config = Config()
    list_exons = {}
    dict_junctions = {}
    list_introns = {}
    dict_of_genes = majiq_io.retrieve_db_genes(majiq_config.outDir)
    for gne_id, gene_obj in dict_of_genes.items():
        list_exons[gne_id] = []
        dict_junctions[gne_id] = {}
        majiq_io.retrieve_db_info(gne_id, majiq_config.outDir, list_exons[gne_id],
                                  dict_junctions[gne_id], None)

    for gne_id, gene_obj in dict_of_genes.items():
        detect_exons(dict_junctions[gne_id], list_exons[gne_id])
        list_introns.update({(gne_id, gene_obj['chromosome'], gene_obj['strand'],
                              ex.end + 1, list_exons[gne_id][id_ex + 1].start - 1): 0
                             for id_ex, ex in enumerate(list_exons[gne_id][:-1])
                             if ex.end + 1 < list_exons[gne_id][id_ex + 1].start - 1})

    for is_junc_file, fname, name in file_list:
        logger.info('READ introns from %s' % fname)
        find_introns(fname, list_introns, majiq_config.min_intronic_cov, process_conf.queue, gname=name)


def parse_denovo_elements(pipe_self, logger):

    majiq_config = Config()
    min_experiments = 2 if majiq_config.min_exp == -1 else majiq_config.min_exp
    dict_of_genes = majiq_io.retrieve_db_genes(majiq_config.outDir)
    list_exons = {}
    dict_junctions = {}
    list_introns = {}

    # nthreads = min(pipe_self.nthreads, len(majiq_config.sam_list))

    lchnksize = math.ceil(len(majiq_config.sam_list) / pipe_self.nthreads)
    nthreads = min(pipe_self.nthreads, math.ceil(len(majiq_config.sam_list) / lchnksize))
    logger.info("Create %s processes" % nthreads)
    pipe_self.lock = [mp.Lock() for _ in range(nthreads)]
    print (nthreads, len(majiq_config.sam_list), lchnksize)
    pool1 = mp.Pool(processes=nthreads, initializer=majiq_multi.process_conf, initargs=[find_new_junctions, pipe_self],
                    maxtasksperchild=1)
    pool2 = mp.Pool(processes=nthreads, initializer=majiq_multi.process_conf, initargs=[find_new_introns, pipe_self],
                    maxtasksperchild=1)

    [xx.acquire() for xx in pipe_self.lock]
    pool1.imap_unordered(majiq_multi.process_wrapper,
                        majiq_multi.chunks(majiq_config.juncfile_list, lchnksize, range(nthreads)))
    pool1.close()


    group_names = {}
    group_lens = []
    for xidx, xx in enumerate(majiq_config.tissue_repl.keys()):
        group_names[xx] = xidx
        group_lens.append(len(xx))
    group_lens = np.array(group_lens)

    with h5py.File(get_build_temp_db_filename(majiq_config.outDir), "r+") as db_f:
        denovo_junctions = {}
        queue_manager(db_f, pipe_self.lock, pipe_self.queue, num_chunks=nthreads, logger=logger,
                      junctions=denovo_junctions, group_names=group_names)
        pool1.join()
        for jj, vv in denovo_junctions.items():
            #vv /= group_lens
            if np.any(vv >= majiq_config.min_denovo):
                #dict_junctions[jj[0]].append(create_junction(jj[1], jj[2], jj[0], 0, annot=False, intron=False))
                majiq_io.dump_junctions(db_f, jj[0], jj[1], jj[2], annot=False)

    [xx.acquire() for xx in pipe_self.lock]
    pool2.imap_unordered(majiq_multi.process_wrapper,
                         majiq_multi.chunks(majiq_config.juncfile_list, lchnksize, range(nthreads)))
    pool2.close()

    with h5py.File(get_build_temp_db_filename(majiq_config.outDir), "r+") as db_f:
        denovo_introns = {}
        queue_manager(db_f, pipe_self.lock, pipe_self.queue, num_chunks=nthreads, logger=logger,
                      introns=denovo_introns, group_names=group_names)
        pool2.join()
        for info, vv in denovo_introns.items():
            if np.any(vv >= min_experiments):
                majiq_io.dump_intron(db_f, info[0], info[1], info[2], annot=False)

    for gne_id, gene_obj in dict_of_genes.items():
        list_exons[gne_id] = []
        dict_junctions[gne_id] = {}
        list_introns[gne_id] = []
        majiq_io.retrieve_db_info(gne_id, majiq_config.outDir, list_exons[gne_id],
                                  dict_junctions[gne_id],
                                  list_introns[gne_id])
        detect_exons(dict_junctions[gne_id], list_exons[gne_id])

    init_splicegraph(get_builder_splicegraph_filename(majiq_config.outDir))

    gene_to_splicegraph(dict_of_genes, dict_junctions, list_exons, list_introns, majiq_config, None)


def parsing_files(sam_file_list, chnk, process_conf, logger):

    majiq_config = Config()

    dict_of_genes = majiq_io.retrieve_db_genes(majiq_config.outDir)
    ngenes = len(dict_of_genes)
    list_exons = {}
    dict_junctions = {}
    list_introns = {}
    effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1

    for gne_id, gene_obj in dict_of_genes.items():
        list_exons[gne_id] = []
        dict_junctions[gne_id] = {}
        list_introns[gne_id] = []
        majiq_io.retrieve_db_info(gne_id, majiq_config.outDir, list_exons[gne_id],
                                  dict_junctions[gne_id],
                                  list_introns[gne_id], default_index=0)
        detect_exons(dict_junctions[gne_id], list_exons[gne_id])

    for sam_file in sam_file_list:
        logger.info("[%s] Starting new file" % sam_file)
        loop_id = sam_file
        samfl = majiq_io_bam.open_rnaseq("%s/%s.bam" % (majiq_config.sam_dir, sam_file))
        junc_mtrx = [[0] * effective_len]

        gc_pairs = {'GC': [], 'COV': []}
        gc_matrx = [] if majiq_config.gcnorm else None

        for gne_idx, (gne_id, gene_obj) in enumerate(dict_of_genes.items()):
            if gne_idx % 50 == 0:
                logger.info("[%s] Progress %s/%s" % (loop_id, gne_idx, ngenes))

            logger.debug("[%s] Reading BAM files" % gne_id)
            gene_reads = majiq_io_bam.read_sam_or_bam(gene_obj, samfl, junc_mtrx, dict_junctions[gne_id],
                                                      list_exons[gne_id], list_introns[gne_id],
                                                      info_msg=loop_id, logging=logger)
            gene_obj['nreads'] = gene_reads
            if gene_reads == 0:
                continue

            if majiq_config.gcnorm:
                for ex in gene_obj.get_exon_list():
                    if ex.get_gc_content() > 0 and ex.get_coverage() > 0:
                        gc_pairs['GC'].append(ex.get_gc_content())
                        gc_pairs['COV'].append(ex.get_coverage())

        majiq_io_bam.close_rnaseq(samfl)

        junc_mtrx = np.array(junc_mtrx)
        update_splicegraph_junctions(dict_junctions, junc_mtrx, majiq_config.outDir, sam_file, process_conf.lock)

        indx = np.arange(junc_mtrx.shape[0])[junc_mtrx.sum(axis=1) >= majiq_config.minreads]

        logger.debug("[%s] Fitting NB function with constitutive events..." % sam_file)
        fitfunc_r = fit_nb(junc_mtrx[indx, :], "%s/nbfit" % majiq_config.outDir, logger=logger)

        with h5py.File('%s/%s.majiq' % (majiq_config.outDir, sam_file), 'w') as out_f:
            #TODO: keep in mem and fix later

            out_f.create_dataset(JUNCTIONS_DATASET_NAME, (5000, majiq_config.m), maxshape=(None, majiq_config.m))
            out_f.create_dataset('junc_cov', (5000, 2), maxshape=(None, 2))

            out_f.attrs['m_samples'] = process_conf.m
            out_f.attrs['sample_id'] = sam_file
            out_f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
            out_f.attrs['VERSION'] = VERSION
            out_f.attrs['lsv_idx'] = 0
            out_f.attrs['num_lsvs'] = 0
            out_f.attrs['genome'] = majiq_config.genome
            out_f.attrs['one_over_r'] = fitfunc_r

            for gne_idx, (gne_id, gene_obj) in enumerate(dict_of_genes.items()):
                if gene_obj['nreads'] == 0:
                    continue

                detect_lsvs(list_exons[gne_id], junc_mtrx, fitfunc_r, gne_id, gene_obj['chromosome'],
                            gene_obj['strand'], majiq_config, out_f)

                for jj in dict_junctions[gne_id].values():
                    jj.index = -1

            n_juncs = out_f.attrs['lsv_idx']
            shp = out_f[JUNCTIONS_DATASET_NAME].shape
            out_f[JUNCTIONS_DATASET_NAME].resize((n_juncs, shp[1]))



        # #TODO: GC CONTENT
        # majiq_norm.mark_stacks(junc_mtrx, fitfunc_r, majiq_config.pvalue_limit)

        # out_f.create_dataset(JUNCTIONS_DATASET_NAME, data=junc_mtrx, compression='gzip', compression_opts=9)
        # out_f.create_dataset(JUNCTIONS_DATASET_NAME, data=junc_mtrx, compression='gzip', compression_opts=9)
        # if majiq_config.gcnorm:
        #     gc_matrx = np.array(gc_matrx)
        #     out_f.create_dataset(JUNCTIONS_GC_CONTENT, data=gc_matrx, compression='gzip', compression_opts=9,
        #                          dtype=np.float)
        #     factor, meanbins = gc_factor_calculation(gc_pairs, nbins=10)
        #     out_f.attrs['gc_values'] = (factor, meanbins)



class Builder(BasicPipeline):

    def run(self):
        if self.simplify is not None and len(self.simplify) not in (0, 2):
            raise RuntimeError('Simplify requires 2 values type of junctions afected and E(PSI) threshold.')
        if not os.path.exists(self.conf):
            raise RuntimeError("Config file %s does not exist" % self.conf)
        majiq_config = Config(self.conf, self)
        self.builder(majiq_config)

    def builder(self, majiq_config):

        logger = majiq_logger.get_logger("%s/majiq.log" % majiq_config.outDir, silent=False,
                                        debug=self.debug)
        logger.info("Majiq Build v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))
        self.queue = mp.Queue()

        if self.prebam:
            p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
                           args=(majiq_io.parse_annot, [self.transcripts, majiq_config.outDir],
                                 '%s/tmp' % majiq_config.outDir, 0))

            logger.info("... waiting gff3 parsing")
            p.start()
            p.join()

            p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
                           args=(parse_denovo_elements, [self], majiq_config.outDir, 0))
            logger.info("... retrieve denovo features")
            p.start()
            p.join()



            #parse_denovo_elements(majiq_config.juncfile_list, 0, None, logger)

        logger.info("Parsing seq files")

        if self.nthreads > 1:
            self.lock = mp.Lock()
            pool = mp.Pool(processes=self.nthreads,
                           initializer=majiq_multi.process_conf,
                           initargs=[parsing_files, self],
                           maxtasksperchild=1)

            lchnksize = max(int(len(majiq_config.sam_list)/self.nthreads), 1) + 1
            pool.imap_unordered(majiq_multi.process_wrapper,
                                majiq_multi.chunks(majiq_config.sam_list, lchnksize, range(self.nthreads)))
            pool.close()
            pool.join()
        else:
            parsing_files(majiq_config.sam_list, 0, process_conf=self, logger=logger)

        for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            with h5py.File('%s/%s.majiq' % (majiq_config.outDir, sam_file), 'r+') as f:
                logger.info('%s LSVs found in %s' % (f.attrs['num_lsvs'], sam_file))

        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
