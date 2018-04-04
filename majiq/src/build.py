import multiprocessing as mp
import os
import sys

import psutil

import majiq.src.io as majiq_io
# from majiq.src.io_bam import find_introns
# import majiq.src.io_bam as majiq_io_bam
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper
from majiq.src.polyfitnb import fit_nb
from voila.api import SpliceGraph
import majiq.src.multiproc as majiq_multi
# from majiq.grimoire.exon import detect_exons, expand_introns
# from majiq.grimoire.lsv import detect_lsvs, sample_junctions
import math
import datetime

from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.config import Config
from majiq.src.constants import *
import majiq.src.logger as majiq_logger
from majiq.src.voila_wrapper import generate_splicegraph, update_splicegraph_junction

from majiq.src.internals.seq_parse import core_build

def build(args):
    pipeline_run(Builder(args))

class Builder(BasicPipeline):

    def store_results(self, output, results, msg_type, extra={}):
        if msg_type == QUEUE_MESSAGE_BUILD_JUNCTION:

            info_junc = results[:-1]
            gidx = extra['group_names'][results[-1]]
            try:
                extra['gen_dict'][info_junc][gidx] += 1
            except KeyError:
                extra['gen_dict'][info_junc] = np.zeros(len(extra['group_names']))
                extra['gen_dict'][info_junc][gidx] = 1

            if extra['gen_dict'][info_junc][gidx] == extra['min_experients'] and info_junc not in extra['found']:
                try:
                    extra['elem_dict'][info_junc[0]].append([info_junc[1], info_junc[2], 0, J_TYPE])
                except KeyError:
                    extra['elem_dict'][info_junc[0]] = [[info_junc[1], info_junc[2], 0, J_TYPE]]
                    extra['found'][info_junc] = 1

        elif msg_type == QUEUE_MESSAGE_BUILD_INTRON:
            info_intron =results[:-1]
            gidx = extra['group_names'][results[-1]]
            try:
                extra['gen_dict'][info_intron][gidx] += 1
            except KeyError:
                extra['gen_dict'][info_intron] = np.zeros(len(extra['group_names']))
                extra['gen_dict'][info_intron][gidx] = 1

            if extra['gen_dict'][info_intron][gidx] == extra['min_experients'] and info_intron not in extra['found']:
                try:
                    extra['elem_dict'][info_intron[0]].append([info_intron[1], info_intron[2], 0, IR_TYPE])
                except KeyError:
                    extra['elem_dict'][info_intron[0]] = [[info_intron[1], info_intron[2], 0, IR_TYPE]]
                    extra['found'][info_intron] = 1

        elif msg_type == QUEUE_MESSAGE_SPLICEGRAPH:
            info = results
            update_splicegraph_junction(output, info[0], info[1], info[2], info[3], info[4])

    def run(self):
        # if self.simplify is not None and len(self.simplify) not in (0, 2):
        #     raise RuntimeError('Simplify requires 2 values type of junctions afected and E(PSI) threshold.')
        if not os.path.exists(self.conf):
            raise RuntimeError("Config file %s does not exist" % self.conf)
        majiq_config = Config(self.conf, self)
        self.builder(majiq_config)

    def builder(self, majiq_config):

        logger = majiq_logger.get_logger("%s/majiq.log" % majiq_config.outDir, silent=False, debug=self.debug)
        logger.info("Majiq Build v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))

        core_build(self.transcripts, majiq_config.sam_list, majiq_config, logger)

        # self.elem_dict =
        # self.genes_dict =
        #
        # mem_allocated = int(psutil.Process().memory_info().rss)/(1024**2)
        # logger.info("PRE DB %.2f MB" % mem_allocated)
        #
        # if self.use_db:
        #     logger.info("Loading previously generated db %s" % self.transcripts)
        #     majiq_io.load_db(self.transcripts, self.elem_dict, self.genes_dict)
        # else:
        #
        #     #majiq_io.parse_annot(self.transcripts, majiq_config.outDir, self.elem_dict, self.genes_dict, logger)
        #     p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
        #                    args=(majiq_io.parse_annot, [self.transcripts, majiq_config.outDir,
        #                                                 self.elem_dict, self.genes_dict],
        #                          '%s/tmp' % majiq_config.outDir, 0))
        #
        #     logger.info("... waiting gff3 parsing")
        #     p.start()
        #     p.join()
        #
        # # p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
        # #                args=(parse_denovo_elements, [self], majiq_config.outDir, 0))
        # # logger.info("... retrieve denovo features")
        # #
        # # p.start()
        # # p.join()
        # mem_allocated = int(psutil.Process().memory_info().rss)/(1024**2)
        # logger.info("PRE parsed %.2f MB" % mem_allocated)
        #
        # self.parse_denovo_elements(logger)
        #
        # mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        # logger.info("POST parsed %.2f MB" % mem_allocated)
        #
        # logger.info("Parsing seq files")
        # # if self.nthreads > 1:
        # #     nthreads = min(self.nthreads, len(majiq_config.sam_list))
        # #     self.lock = [mp.Lock() for _ in range(nthreads)]
        # #     [xx.acquire() for xx in self.lock]
        # #     pool = mp.Pool(processes=nthreads,
        # #                    initializer=majiq_multi.process_conf,
        # #                    initargs=[parsing_files, self],
        # #                    maxtasksperchild=1)
        # #
        # #     pool.imap_unordered(majiq_multi.process_wrapper,
        # #                         majiq_multi.chunks(majiq_config.sam_list, nthreads))
        # #     pool.close()
        # #     generate_splicegraph(majiq_config, self.elem_dict, self.genes_dict)
        # #     with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'r+') as sg:
        # #        queue_manager(sg, self.lock, self.queue, num_chunks=nthreads, func=self.store_results, logger=logger)
        # #
        # #     pool.join()
        # # else:
        # #     parsing_files(majiq_config.sam_list, 0, conf=self, logger=logger)
        #
        # # for exp_idx, sam_file in enumerate(majiq_config.sam_list):
        #     # with h5py.File('%s/%s.majiq' % (majiq_config.outDir, sam_file), 'r+') as f:
        #     #     logger.info('%s LSVs found in %s' % (f.attrs['num_lsvs'], sam_file))

        if self.mem_profile:
            mem_allocated = int(psutil.Process().memory_info().rss)/(1024**2)
            logger.info("Max Memory used %.2f MB" % mem_allocated)
        logger.info("MAJIQ Builder is ended succesfully!")

