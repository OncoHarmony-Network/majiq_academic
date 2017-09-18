import multiprocessing as mp
import os
import sys

import h5py

import majiq.src.io as majiq_io
import majiq.src.io_bam as io_bam
import majiq.src.logger as majiq_logger
import majiq.src.multiproc as majiq_multi
import majiq.src.normalize as majiq_norm
from majiq.grimoire.exon import detect_exons
from majiq.grimoire.lsv import detect_lsvs
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.config import Config
from majiq.src.constants import *
from majiq.src.polyfitnb import fit_nb


def build(args):
    pipeline_run(Builder(args))


def merging_files(list_of_genes, chnk, process_conf, logger):
    majiq_config = Config()
    total = len(list_of_genes)
    logger.debug("[%s] List of genes" % ",".join(list_of_genes))
    for gne_idx, gne_id in enumerate(list_of_genes):
        list_exons = []
        dict_junctions = {}
        fitfunc_r = []
        if gne_idx % 50 == 0:
            logger.info("[%s] Progress %s/%s" % (chnk, gne_idx, total))

        loop_id = '%s - %s' % (chnk, gne_id)
        logger.debug("[%s] Retrieving gene" % loop_id)

        """Retrieve annot junctions from exon """
        majiq_io.retrieve_db_info(gne_id, majiq_config.outDir, list_exons, dict_junctions)

        """Retrieve quantifications and denovo juncs"""
        junc_mtrx = majiq_io.get_covered_junctions(gne_id, dict_junctions, list_exons, fitfunc_r, majiq_config.sam_list,
                                                   majiq_config.readLen, majiq_config.outDir)

        """detect exons"""
        detect_exons(dict_junctions, list_exons)

        """IR adaptation"""

        if len(list_exons) > 1 and majiq_config.ir:
            if majiq_config.ir:
                logger.debug("[%s] Detecting intron retention events" % gne_id)
                io_bam.rnaseq_intron_retention(gene_obj, list_exons[gne_id], majiq_config.sam_list, junc_mtrx,
                                               out_junctions=junc_list, logging=logger)

        detect_lsvs(list_exons, junc_mtrx, fitfunc_r, process_conf.lock, gne_id, '+', majiq_config)

        del list_exons


def parsing_files(sam_file_list, chnk, process_conf, logger):

    majiq_config = Config()
    dict_of_genes = majiq_io.retrieve_db_genes(majiq_config.outDir)
    ngenes = len(dict_of_genes)
    if majiq_config.ir:
        list_exons = {}
        for gne_id, gene_obj in dict_of_genes.items():
            list_exons[gne_id] = []
            majiq_io.retrieve_db_info(gne_id, majiq_config.outDir, list_exons[gne_id], None)

    for sam_file in sam_file_list:
        logger.info("[%s] Starting new file" % sam_file)
        loop_id = sam_file
        out_f = h5py.File(get_builder_temp_majiq_filename(majiq_config.outDir, sam_file), 'w')
        samfl = io_bam.open_rnaseq("%s/%s.bam" % (majiq_config.sam_dir, sam_file))
        gc_pairs = {'GC': [], 'COV': []}

        junc_mtrx = []
        gc_matrx = [] if majiq_config.gcnorm else None

        for gne_idx, (gne_id, gene_obj) in enumerate(dict_of_genes.items()):
            if gne_idx % 50 == 0:
                logger.info("[%s] Progress %s/%s" % (loop_id, gne_idx, ngenes))

            logger.debug("[%s] Reading BAM files" % gne_id)
            junc_list = []
            gene_reads = io_bam.read_sam_or_bam(gene_obj, samfl, junc_mtrx, out_junctions=junc_list,
                                                info_msg=loop_id, logging=logger)
            if gene_reads == 0:
                continue

            if majiq_config.gcnorm:
                for ex in gene_obj.get_exon_list():
                    if ex.get_gc_content() > 0 and ex.get_coverage() > 0:
                        gc_pairs['GC'].append(ex.get_gc_content())
                        gc_pairs['COV'].append(ex.get_coverage())

            for jnc in junc_list:
                majiq_io.junction_to_tmp(gne_id, jnc, out_f)

        io_bam.close_rnaseq(samfl)
        junc_mtrx = np.array(junc_mtrx)
        indx = np.arange(junc_mtrx.shape[0])[junc_mtrx.sum(axis=1) >= majiq_config.minreads]
        logger.debug("[%s] Fitting NB function with constitutive events..." % sam_file)

        fitfunc_r = fit_nb(junc_mtrx[indx, :], "%s/nbfit" % majiq_config.outDir, logger=logger)

        #TODO: GC CONTENT
        majiq_norm.mark_stacks(junc_mtrx, fitfunc_r, majiq_config.pvalue_limit)
        out_f.attrs['one_over_r'] = fitfunc_r

        out_f.create_dataset(JUNCTIONS_DATASET_NAME, data=junc_mtrx, compression='gzip', compression_opts=9)
        # out_f.create_dataset(JUNCTIONS_DATASET_NAME, data=junc_mtrx, compression='gzip', compression_opts=9)
        # if majiq_config.gcnorm:
        #     gc_matrx = np.array(gc_matrx)
        #     out_f.create_dataset(JUNCTIONS_GC_CONTENT, data=gc_matrx, compression='gzip', compression_opts=9,
        #                          dtype=np.float)
        #     factor, meanbins = gc_factor_calculation(gc_pairs, nbins=10)
        #     out_f.attrs['gc_values'] = (factor, meanbins)

        out_f.close()

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

        if self.prebam:

            p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
                           args=(majiq_io.parse_annot, [self.transcripts, majiq_config.outDir, majiq_config.juncfile_list,
                                                        majiq_config.denovo, majiq_config.min_denovo],
                                 '%s/tmp' % majiq_config.outDir, 'db', 0))
            logger.info("... waiting gff3 parsing")
            p.start()
            p.join()

            if self.nthreads > 1:
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

        list_of_genes = majiq_io.get_list_of_genes(majiq_config.outDir)
        lchnksize = max(len(list_of_genes)/self.nthreads, 1) + 1

        self.lock = [mp.Lock() for _ in range(len(majiq_config.sam_list))]
        self.lock.append(mp.Lock())

        pool = mp.Pool(processes=self.nthreads,
                       initializer=majiq_multi.process_conf,
                       initargs=[merging_files, self],
                       maxtasksperchild=1)

        #init_splicegraph(get_builder_splicegraph_filename(majiq_config.outDir))

        for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            majiq_io.init_majiq_file(sam_file, majiq_config.outDir, genome=majiq_config.genome, msamples=majiq_config.m)

        pool.imap_unordered(majiq_multi.process_wrapper,
                            majiq_multi.chunks(list_of_genes, lchnksize, range(self.nthreads)))
        pool.close()
        pool.join()

        for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            with h5py.File('%s/%s.majiq' % (majiq_config.outDir, sam_file), 'r+') as f:
                n_juncs = f.attrs['lsv_idx']
                shp = f[JUNCTIONS_DATASET_NAME].shape
                f[JUNCTIONS_DATASET_NAME].resize((n_juncs, shp[1]))
                logger.info('%s LSVs found in %s' % (f.attrs['num_lsvs'], sam_file))

        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
