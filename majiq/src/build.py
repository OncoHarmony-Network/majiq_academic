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
from majiq.src.sample import create_lt

from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.config import Config
from majiq.src.constants import *
from majiq.src.polyfitnb import fit_nb
from majiq.src.voila_wrapper import generate_splicegraph
from voila.api import SpliceGraph

from quicksect import IntervalNode, Interval, IntervalTree

import math
import datetime


def build(args):
    pipeline_run(Builder(args))

def find_new_junctions(file_list, chunk, conf, logger):

    majiq_config = Config()
    list_exons = {}
    dict_junctions = {}
    logger.info('Reading DB')

    for gne_id, gene_obj in conf.genes_dict.items():
        dict_junctions[gne_id] = {}
        list_exons[gne_id] = []

        majiq_io.from_matrix_to_objects(gne_id, conf.elem_dict[gne_id], dict_junctions[gne_id], list_exons[gne_id])
        detect_exons(dict_junctions[gne_id], list_exons[gne_id])

    td = create_lt(conf.genes_dict)
    for is_junc_file, fname, name in file_list:
        logger.info('READ JUNCS from %s, %s' % (fname, majiq_config.strand_specific))
        # read_juncs2(fname, is_junc_file, list_exons, conf.genes_dict, dict_junctions,
        #            majiq_config.strand_specific, process_conf.queue, gname=name)

        read_juncs(fname, is_junc_file, list_exons, conf.genes_dict, td, dict_junctions,
                   majiq_config.strand_specific, conf.queue, gname=name)



def find_new_introns(file_list, chunk, conf, logger):

    majiq_config = Config()
    list_introns = {}

    num_bins = 10
    logger.info('Reading DB')

    for gne_id, gene_obj in conf.genes_dict.items():
        list_exons = []
        dict_junctions = {}
        introns = []
        majiq_io.from_matrix_to_objects(gne_id, conf.elem_dict[gne_id], dict_junctions, list_exons, introns)

        detect_exons(dict_junctions, list_exons)

        range_introns = IntervalTree()
        nir = len(introns)
        [range_introns.add(xx.start, xx.end+1) for xx in introns]
        del introns

        list_introns = dict()
        for id_ex, ex in enumerate(list_exons[:-1]):
            if (ex.end + 3 < list_exons[id_ex + 1].start - 1 and ex.end != -1
               and list_exons[id_ex + 1].start != -1) :

                if nir > 0 and len(range_introns.search(ex.end + 1, ex.end + 2)) != 0:
                    continue

                ilen = (list_exons[id_ex + 1].start - 1) - (ex.end + 1)
                nchunks = 1 if ilen <= MIN_INTRON_LEN else num_bins
                chunk_len = int(ilen / nchunks)+1

                try:
                    list_introns[gene_obj['chromosome']].add(ex.end + MIN_BP_OVERLAP + 1,
                                                             list_exons[id_ex + 1].start - (MIN_BP_OVERLAP + 1),
                                                             (gne_id, gene_obj['strand'], nchunks, chunk_len))

                except KeyError:
                    list_introns[gene_obj['chromosome']] = IntervalTree()
                    list_introns[gene_obj['chromosome']].add(ex.end + MIN_BP_OVERLAP + 1,
                                                             list_exons[id_ex + 1].start - (MIN_BP_OVERLAP + 1),
                                                             (gne_id, gene_obj['strand'], nchunks, chunk_len))

        del range_introns

    for is_junc_file, fname, name in file_list:
        logger.info('READ introns from %s' % fname)
        find_introns(fname, list_introns, majiq_config.min_intronic_cov, conf.queue, gname=name)


def parse_denovo_elements(pself, logger):

    majiq_config = Config()
    min_experiments = 2 if majiq_config.min_exp == -1 else majiq_config.min_exp

    elem_dict = {}

    if majiq_config.denovo:
        nthreads = min(pself.nthreads, len(majiq_config.sam_list))
        logger.info("Create %s processes" % nthreads)
        pself.lock = [mp.Lock() for _ in range(nthreads)]

        pool1 = mp.Pool(processes=nthreads, initializer=majiq_multi.process_conf,
                        initargs=[find_new_junctions, pself], maxtasksperchild=1)
        if majiq_config.ir:
            pool2 = mp.Pool(processes=nthreads, initializer=majiq_multi.process_conf,
                            initargs=[find_new_introns, pself], maxtasksperchild=1)

        [xx.acquire() for xx in pself.lock]
        pool1.imap_unordered(majiq_multi.process_wrapper,
                             majiq_multi.chunks(majiq_config.juncfile_list, nthreads))
        pool1.close()

        group_names = {xx: xidx for xidx, xx in enumerate(majiq_config.tissue_repl.keys())}
        queue_manager(None, pself.lock, pself.queue, num_chunks=nthreads, logger=logger,
                      elem_dict=elem_dict, group_names=group_names, min_experients=min_experiments)
        pool1.join()

        logger.info('Updating DB')
        majiq_io.add_elements_mtrx(elem_dict, pself.elem_dict)

        if majiq_config.ir:
            elem_dict = {}
            [xx.acquire() for xx in pself.lock]
            pool2.imap_unordered(majiq_multi.process_wrapper,
                                 majiq_multi.chunks(majiq_config.juncfile_list, nthreads))
            pool2.close()
            queue_manager(None, pself.lock, pself.queue, num_chunks=nthreads, logger=logger,
                          elem_dict=elem_dict, group_names=group_names, min_experients=min_experiments)
            pool2.join()
            majiq_io.add_elements_mtrx(elem_dict, pself.elem_dict)

    majiq_io.dump_elements(pself.genes_dict, pself.elem_dict, pself.outDir)


def parsing_files(sam_file_list, chnk, conf, logger):

    majiq_config = Config()
    effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1

    list_exons = {}
    dict_junctions = {}
    list_introns = {}
    logger.info('Reading DB')

    ngenes = len(conf.genes_dict)
    loop_step = max(1, int(ngenes/10))

    for gne_id, gene_obj in conf.genes_dict.items():
        dict_junctions[gne_id] = {}
        list_exons[gne_id] = []
        list_introns[gne_id] = []
        majiq_io.from_matrix_to_objects(gne_id, conf.elem_dict[gne_id], dict_junctions[gne_id], list_exons[gne_id],
                                        list_introns[gne_id], default_index=0)
        detect_exons(dict_junctions[gne_id], list_exons[gne_id])
        if majiq_config.ir:
            expand_introns(gne_id, list_introns[gne_id], list_exons[gne_id], dict_junctions[gne_id], default_index=0)

    for sam_file in sam_file_list:
        logger.info("[%s] Starting new file" % sam_file)
        loop_id = sam_file
        samfl = majiq_io_bam.open_rnaseq("%s/%s.bam" % (majiq_config.sam_dir, sam_file))
        junc_mtrx = [[0] * effective_len]

        for gne_idx, (gne_id, gene_obj) in enumerate(conf.genes_dict.items()):
            if gne_idx % loop_step == 0:
                logger.info("[%s] Progress %s/%s" % (loop_id, gne_idx, ngenes))

            logger.debug("[%s] Reading BAM files" % gne_id)
            gene_reads = majiq_io_bam.read_sam_or_bam(gene_obj, samfl, junc_mtrx, dict_junctions[gne_id],
                                                      list_exons[gne_id], list_introns[gne_id],
                                                      info_msg=loop_id, logging=logger)
            # gene_obj['nreads'] = gene_reads
            if gene_reads == 0:
                continue

            for jj in dict_junctions[gne_id].values():
                qm = QueueMessage(QUEUE_MESSAGE_SPLICEGRAPH,
                                  (gne_id, jj.start, jj.end, np.sum(junc_mtrx[jj.index]), sam_file), chnk)
                conf.queue.put(qm, block=True)

        majiq_io_bam.close_rnaseq(samfl)
        junc_mtrx = np.array(junc_mtrx, dtype=np.float)

        indx = np.arange(junc_mtrx.shape[0])[junc_mtrx.sum(axis=1) >= majiq_config.minreads]
        logger.debug("[%s] Fitting NB function with constitutive events..." % sam_file)
        fitfunc_r = fit_nb(junc_mtrx[indx, :], "%s/nbfit" % majiq_config.outDir, logger=logger)

        # with h5py.File('%s/%s.majiq' % (majiq_config.outDir, sam_file), 'w') as out_f:
        #     out_f.attrs['m_samples'] = process_conf.m
        #     out_f.attrs['sample_id'] = sam_file
        #     out_f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
        #     out_f.attrs['VERSION'] = VERSION

        logger.info('Detecting lsvs')
        lsv_type_list = []
        lsv_dict = {}

        detect_lsvs(conf.genes_dict, dict_junctions, list_exons, junc_mtrx, fitfunc_r, majiq_config, lsv_dict,
                    lsv_type_list, logger)

        logger.info('dump samples')
        # vals = sample_junctions(lsv_dict, fitfunc_r, majiq_config)
        fname = '%s/%s.majiq' % (majiq_config.outDir, sam_file)
        majiq_io.dump_lsv_coverage(fname, lsv_dict, lsv_type_list)

        del lsv_type_list
        del lsv_dict


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


        manager = mp.Manager()
        self.elem_dict = manager.dict()
        self.genes_dict = manager.dict()

        if self.use_db:
            logger.info("Loading previously generated db %s" % self.transcripts)
            majiq_io.load_db(self.transcripts, self.elem_dict, self.genes_dict)
        else:
            p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
                           args=(majiq_io.parse_annot, [self.transcripts, majiq_config.outDir,
                                                        self.elem_dict, self.genes_dict],
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
            generate_splicegraph(majiq_config, self.elem_dict, self.genes_dict)
            with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'r+') as sg:
                queue_manager(sg, self.lock, self.queue, num_chunks=nthreads, logger=logger)

            pool.join()
        else:
            parsing_files(majiq_config.sam_list, 0, conf=self, logger=logger)

        # for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            # with h5py.File('%s/%s.majiq' % (majiq_config.outDir, sam_file), 'r+') as f:
            #     logger.info('%s LSVs found in %s' % (f.attrs['num_lsvs'], sam_file))

        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
