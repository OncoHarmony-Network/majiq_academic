#!/usr/bin/python

import argparse
import os
import sys
import traceback
from multiprocessing import Pool, current_process, Process
from majiq.grimoire.gene import recreate_gene_tlb

import majiq.src.analize as analize
import majiq.src.io as majiq_io
from majiq.src.normalize import prepare_gc_content
import majiq.src.utils.utils as majiq_utils
import majiq.src.config as mglobals
import majiq.grimoire.lsv as majiq_lsv

try:
    import cPickle as pickle
except Exception:
    import pickle


def majiq_builder(samfiles_list, chnk, pcr_validation=None, gff_output=None, create_tlb=True, only_rna=False,
                  nondenovo=False, logging=None):

    majiq_utils.monitor('CHUNK %s: INITIAL' % chnk)

    if not logging is None:
        logging.info("Building chunk %s" % chnk)

    temp_dir = "%s/tmp/chunk_%s" % (mglobals.outDir, chnk)
    annot_file = '%s/annot_genes.pkl' % temp_dir
    if not os.path.exists(annot_file):
        return
    temp_file = open(annot_file, 'rb')
    gene_list = pickle.load(temp_file)

    if create_tlb:
        if not logging is None:
            logging.info("[%s] Recreatin Gene TLB" % chnk)
        recreate_gene_tlb(gene_list)

    if not logging is None:
        logging.info("[%s] Reading BAM files" % chnk)
    majiq_utils.monitor('CHUNK %s: PRE SAM' % chnk)
    majiq_io.read_sam_or_bam(samfiles_list, gene_list, chnk,
                             nondenovo=nondenovo, logging=logging)
    if not logging is None:
        logging.info("[%s] Detecting intron retention events" % chnk)
    majiq_io.rnaseq_intron_retention(samfiles_list, gene_list, chnk,
                                     permissive=mglobals.permissive_ir, nondenovo=nondenovo, logging=logging)
    if not logging is None:
        logging.info("[%s] Detecting LSV" % chnk)
    majiq_utils.monitor('CHUNK %s: PRE DETECTION' % chnk)
    lsv, const = analize.lsv_detection(gene_list, chnk, only_real_data=only_rna, logging=logging)

    prepare_gc_content(gene_list, temp_dir)

    if pcr_validation:
        majiq_utils.get_validated_pcr_lsv(lsv, temp_dir)
    if gff_output:
        majiq_lsv.extract_gff(lsv, temp_dir)
    majiq_utils.generate_visualization_output(gene_list, temp_dir)
    if not logging is None:
        logging.info("[%s] Preparing output" % chnk)

    majiq_utils.monitor('CHUNK %s: PRE TABLE' % chnk)
    majiq_utils.prepare_lsv_table(lsv, const, temp_dir)
    majiq_utils.monitor('CHUNK %s: END' % chnk)

    #ANALYZE_DENOVO
    # utils.analyze_denovo_junctions(gene_list, "%s/denovo.pkl" % temp_dir)
    # utils.histogram_for_exon_analysis(gene_list, "%s/ex_lengths.pkl" % temp_dir)


def __parallel_lsv_quant(samfiles_list, chnk, pcr_validation=False, gff_output=None, only_rna=False,
                         nondenovo=False, silent=False, debug=0):
    try:
        print "START child,", current_process().name
        tlogger = majiq_utils.get_logger("%s/%s.majiq.log" % (mglobals.outDir, chnk), silent=silent, debug=debug)
        majiq_builder(samfiles_list, chnk, pcr_validation=pcr_validation,
                      gff_output=gff_output, only_rna=only_rna, nondenovo=nondenovo,
                      logging=tlogger)
        print "END child, ", current_process().name
    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


def __parallel_gff3(transcripts, pcr_filename, nthreads, silent=False, debug=0):

    try:
        print "START child,", current_process().name
        tlogger = majiq_utils.get_logger("%s/db.majiq.log" % mglobals.outDir, silent=silent, debug=debug)
        majiq_io.read_gff(transcripts, pcr_filename, nthreads, logging=tlogger)
        print "END child, ", current_process().name
    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


#########
# MAIN  #
#########


def main(params):

    if not os.path.exists(params.conf):
        raise RuntimeError("Config file %s does not exist" % params.conf)
    mglobals.global_conf_ini(params.conf, params)

    logger = majiq_utils.get_logger("%s/majiq.log" % mglobals.outDir, silent=params.silent, debug=params.debug)
    logger.info("")
    logger.info("Command: %s" % params)

    if not params.onlygather:
        p = Process(target=__parallel_gff3, args=(params.transcripts, params.pcr_filename, params.nthreads))
        logger.info("... waiting gff3 parsing")
        p.start()
        p.join()

        logger.info("Get samfiles")
        sam_list = []
        for exp_idx, exp in enumerate(mglobals.exp_list):
            samfile = "%s/%s.bam" % (mglobals.sam_dir, exp)
            if not os.path.exists(samfile):
                raise RuntimeError("Skipping %s.... not found" % samfile)
                #logger.info("Skipping %s.... not found" % samfile)
                #continue
            baifile = "%s/%s.bam.bai" % (mglobals.sam_dir, exp)
            if not os.path.exists(baifile):
                raise RuntimeError("Skipping %s.... not found ( index file for bam file is required)" % baifile)
                #logger.info("Skipping %s.... not found ( index file for bam file is required)" % baifile)
                #continue
            sam_list.append(samfile)
            mglobals.exp_list[exp_idx] = os.path.split(exp)[1]
        #majiq_io.count_mapped_reads(samfile, exp_idx)
        if len(sam_list) == 0:
            return

        if params.nthreads > 1:
            pool = Pool(processes=params.nthreads, maxtasksperchild=1)
        #logger.info("Scatter in Chromosomes")
        # for chrom in chr_list:
        majiq_utils.monitor('MAIN: PRE RNA')
        for chnk in range(mglobals.num_final_chunks):
            temp_dir = "%s/tmp/chunk_%s" % (mglobals.outDir, chnk)
            majiq_utils.create_if_not_exists(temp_dir)
            if params.nthreads == 1:
                majiq_builder(sam_list, chnk, pcr_validation=params.pcr_filename, gff_output=params.gff_output,
                              only_rna=params.only_rna, nondenovo=params.non_denovo, logging=logger)
            else:

                pool.apply_async(__parallel_lsv_quant, [sam_list, chnk,
                                                        params.pcr_filename,
                                                        params.gff_output,
                                                        params.only_rna,
                                                        params.non_denovo,
                                                        params.silent,
                                                        params.debug])


        if params.nthreads > 1:
            logger.info("... waiting childs")
            pool.close()
            pool.join()

    majiq_utils.monitor('MAIN: PRE GATHER')
    #GATHER
    majiq_utils.gather_files(mglobals.outDir, '', params.gff_output, params.pcr_filename,
                             nthreads=params.nthreads, logger=logger)

    majiq_utils.monitor('MAIN: END')

    mglobals.print_numbers()
    logger.info("End of execution")
