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
from majiq.grimoire.exon import detect_exons, expand_introns
from majiq.grimoire.lsv import detect_lsvs, sample_junctions
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.config import Config
from majiq.src.constants import *
from majiq.src.polyfitnb import fit_nb
from majiq.src.voila_wrapper import generate_splicegraph
from voila.api import SpliceGraph

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
        majiq_io.retrieve_db(gne_id, majiq_config.outDir, dict_junctions[gne_id], list_exons[gne_id], None)
        detect_exons(dict_junctions[gne_id], list_exons[gne_id])

    for is_junc_file, fname, name in file_list:
        logger.info('READ JUNCS from %s, %s' % (fname, majiq_config.strand_specific))
        read_juncs(fname, is_junc_file, list_exons, dict_of_genes, dict_junctions,
                   majiq_config.strand_specific, process_conf.queue, gname=name)

def find_new_introns(file_list, chunk, process_conf, logger):

    majiq_config = Config()
    list_introns = {}
    dict_of_genes = majiq_io.retrieve_db_genes(majiq_config.outDir)
    num_bins = 10
    for gne_id, gene_obj in dict_of_genes.items():
        list_exons = []
        dict_junctions = {}
        introns = []
        majiq_io.retrieve_db(gne_id, majiq_config.outDir, dict_junctions, list_exons, introns)
        range_introns = [range(xx.start, xx.end+1) for xx in introns]
        del introns

        detect_exons(dict_junctions, list_exons)

        for id_ex, ex in enumerate(list_exons[:-1]):
            if (ex.end + 3 < list_exons[id_ex + 1].start - 1 and ex.end != -1 and list_exons[id_ex + 1].start != -1
               and not np.any([(ex.end + 1) in xx for xx in range_introns])):

                ilen = (list_exons[id_ex + 1].start - 1) - (ex.end + 1)
                nchunks = 1 if ilen <= MIN_INTRON_LEN else num_bins
                chunk_len = int(ilen / nchunks)+1

                try:
                    list_introns[gene_obj['chromosome']].append((gne_id, gene_obj['strand'], ex.end + 1,
                                                                list_exons[id_ex + 1].start - 1, nchunks, chunk_len))
                except KeyError:
                    list_introns[gene_obj['chromosome']] = [(gne_id, gene_obj['strand'], ex.end + 1,
                                                                list_exons[id_ex + 1].start - 1, nchunks, chunk_len)]


    for is_junc_file, fname, name in file_list:
        logger.info('READ introns from %s' % fname)
        find_introns(fname, list_introns, majiq_config.min_intronic_cov, process_conf.queue, gname=name)


def parse_denovo_elements(pipe_self, logger):

    majiq_config = Config()
    min_experiments = 2 if majiq_config.min_exp == -1 else majiq_config.min_exp

    elem_dict = {}

    if majiq_config.denovo:
        nthreads = min(pipe_self.nthreads, len(majiq_config.sam_list))
        logger.info("Create %s processes" % nthreads)
        pipe_self.lock = [mp.Lock() for _ in range(nthreads)]

        pool1 = mp.Pool(processes=nthreads, initializer=majiq_multi.process_conf,
                        initargs=[find_new_junctions, pipe_self], maxtasksperchild=1)
        if majiq_config.ir:
            pool2 = mp.Pool(processes=nthreads, initializer=majiq_multi.process_conf,
                            initargs=[find_new_introns, pipe_self], maxtasksperchild=1)

        [xx.acquire() for xx in pipe_self.lock]
        pool1.imap_unordered(majiq_multi.process_wrapper,
                             majiq_multi.chunks(majiq_config.juncfile_list, nthreads))
        pool1.close()
        group_names = {xx: xidx for xidx, xx in enumerate(majiq_config.tissue_repl.keys())}
        queue_manager(None, pipe_self.lock, pipe_self.queue, num_chunks=nthreads, logger=logger,
                      elem_dict=elem_dict, group_names=group_names, min_experients=min_experiments)
        pool1.join()

        if majiq_config.ir:
            [xx.acquire() for xx in pipe_self.lock]
            pool2.imap_unordered(majiq_multi.process_wrapper,
                                 majiq_multi.chunks(majiq_config.juncfile_list, nthreads))
            pool2.close()
            queue_manager(None, pipe_self.lock, pipe_self.queue, num_chunks=nthreads, logger=logger,
                          elem_dict=elem_dict, group_names=group_names)
            pool2.join()


def parsing_files(sam_file_list, chnk, process_conf, logger):

    majiq_config = Config()

    # dict_of_genes = majiq_io.retrieve_db_genes(majiq_config.outDir)

    effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1

    list_exons = []
    dict_junctions = {}
    list_introns = []
    dict_of_genes = majiq_io.retrieve(majiq_config.outDir, dict_junctions, list_exons,
                                      list_introns, default_index=0)
    ngenes = len(dict_of_genes)

    for gne_id, gene_obj in dict_of_genes.items():

        detect_exons(dict_junctions[gne_id], list_exons[gne_id])
        if majiq_config.ir:
            expand_introns(gne_id, list_introns[gne_id], list_exons[gne_id], dict_junctions[gne_id], default_index=0)
    #     list_exons[gne_id] = []
    #     dict_junctions[gne_id] = {}
    #     list_introns[gne_id] = []
    #     majiq_io.retrieve_db(gne_id, majiq_config.outDir, dict_junctions[gne_id], list_exons[gne_id],
    #                               list_introns[gne_id], default_index=0)
    #     detect_exons(dict_junctions[gne_id], list_exons[gne_id])
    #     if majiq_config.ir:
    #         expand_introns(gne_id, list_introns[gne_id], list_exons[gne_id], dict_junctions[gne_id], default_index=0)

    for sam_file in sam_file_list:
        logger.info("[%s] Starting new file" % sam_file)
        loop_id = sam_file
        samfl = majiq_io_bam.open_rnaseq("%s/%s.bam" % (majiq_config.sam_dir, sam_file))
        junc_mtrx = [[0] * effective_len]

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
            for jj in dict_junctions[gne_id].values():
                qm = QueueMessage(QUEUE_MESSAGE_SPLICEGRAPH,
                                  (gne_id, jj.start, jj.end, np.sum(junc_mtrx[jj.index]), sam_file), chnk)
                process_conf.queue.put(qm, block=True)

        majiq_io_bam.close_rnaseq(samfl)
        junc_mtrx = np.array(junc_mtrx, dtype=np.float)

        indx = np.arange(junc_mtrx.shape[0])[junc_mtrx.sum(axis=1) >= majiq_config.minreads]
        logger.debug("[%s] Fitting NB function with constitutive events..." % sam_file)
        fitfunc_r = fit_nb(junc_mtrx[indx, :], "%s/nbfit" % majiq_config.outDir, logger=logger)

        with h5py.File('%s/%s.majiq' % (majiq_config.outDir, sam_file), 'w') as out_f:
            out_f.attrs['m_samples'] = process_conf.m
            out_f.attrs['sample_id'] = sam_file
            out_f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
            out_f.attrs['VERSION'] = VERSION
            out_f.attrs['lsv_idx'] = 0
            out_f.attrs['num_lsvs'] = 0
            out_f.attrs['genome'] = majiq_config.genome
            out_f.attrs['one_over_r'] = fitfunc_r

            logger.info('Detecting lsvs')
            np_jjlist = [np.zeros(effective_len)]

            detect_lsvs(dict_of_genes, dict_junctions, list_exons, junc_mtrx, fitfunc_r, majiq_config, out_f,
                        np_jjlist, logger)

            logger.info('dump samples')
            vals = sample_junctions(np.array(np_jjlist), fitfunc_r, majiq_config)
            majiq_io.dump_lsv_coverage(out_f, vals[0], vals[1])
            del np_jjlist


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

        if not self.use_db:
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

        logger.info("Parsing seq files")
        if self.nthreads > 1:
            nthreads = min(self.nthreads, len(majiq_config.sam_list))
            self.lock = [mp.Lock() for _ in range(nthreads)]
            [xx.acquire() for xx in self.lock]
            pool = mp.Pool(processes=nthreads,
                           initializer=majiq_multi.process_conf,
                           initargs=[parsing_files, self],
                           maxtasksperchild=1)

            pool.imap_unordered(majiq_multi.process_wrapper,
                                majiq_multi.chunks(majiq_config.sam_list, nthreads))
            pool.close()
            generate_splicegraph(majiq_config)
            with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'r+') as sg:
                queue_manager(sg, self.lock, self.queue, num_chunks=nthreads, logger=logger)

            pool.join()
        else:
            parsing_files(majiq_config.sam_list, 0, process_conf=self, logger=logger)

        for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            with h5py.File('%s/%s.majiq' % (majiq_config.outDir, sam_file), 'r+') as f:
                logger.info('%s LSVs found in %s' % (f.attrs['num_lsvs'], sam_file))

        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
