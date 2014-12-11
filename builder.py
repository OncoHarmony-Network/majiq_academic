#!/usr/bin/python

#standard python libraries

import argparse
import os
import sys
import traceback
from multiprocessing import Pool, current_process, Process
import grimoire.analize as analize
import grimoire.rnaseq_io as majiq_io
import grimoire.utils.utils as utils
import grimoire.mglobals as mglobals
import grimoire.lsv as majiq_lsv

try:
    import cPickle as pickle
except Exception:
    import pickle


def majiq_builder(samfiles_list, chrom, pcr_validation=None, gff_output=None, create_tlb=True, only_rna=False,
                  nondenovo=False, logging=None):

    if not logging is None:
        logging.info("Building for chromosome %s" % chrom)

    temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
    temp_file = open('%s/annot_genes.pkl' % temp_dir, 'rb')
    gene_list = pickle.load(temp_file)

    if create_tlb:
        if not logging is None:
            logging.info("[%s] Recreatin Gene TLB" % chrom)
        utils.recreate_gene_tlb(gene_list)

    if not logging is None:
        logging.info("[%s] Reading BAM files" % chrom)
    majiq_io.read_sam_or_bam(samfiles_list, gene_list, mglobals.readLen, chrom,
                             nondenovo=nondenovo, logging=logging)
    if not logging is None:
        logging.info("[%s] Detecting LSV" % chrom)
    lsv, const = analize.lsv_detection(gene_list, chrom, only_real_data=only_rna, logging=logging)

    utils.prepare_gc_content(gene_list, temp_dir)

    if pcr_validation:
        utils.get_validated_pcr_lsv(lsv, temp_dir)
    if gff_output:
        majiq_lsv.extract_gff(lsv, temp_dir)
    utils.generate_visualization_output(gene_list, temp_dir)
    if not logging is None:
        logging.info("[%s] Preparing output" % chrom)

    utils.prepare_lsv_table(lsv, const, temp_dir)


def __parallel_lsv_quant(samfiles_list, chrom, pcr_validation=False, gff_output=None, silent=False, debug=0):

    try:
        print "START child,", current_process().name
        tlogger = utils.get_logger("%s/%s.majiq.log" % (mglobals.outDir, chrom), silent=silent, debug=debug)
        majiq_builder(samfiles_list, chrom, pcr_validation, gff_output, create_tlb=True, logging=tlogger)
        print "END child, ", current_process().name
    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


def __parallel_gff3(transcripts, pcr_filename, output, silent=False, debug=0):

    try:
        print "START child,", current_process().name
        tlogger = utils.get_logger("%s/db.majiq.log" % mglobals.outDir, silent=silent, debug=debug)
        majiq_io.read_gff(transcripts, pcr_filename, logging=tlogger)
        print "END child, ", current_process().name
    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


def _new_subparser():
    return argparse.ArgumentParser(add_help=False)


def _generate_parser():
    dir_or_paths = _new_subparser()

    parser = argparse.ArgumentParser(parents=[dir_or_paths])
    parser.add_argument('transcripts', action="store", help='read file in SAM format')
    parser.add_argument('-conf', default=None, required=True, help='Provide study configuration file with all '
                                                                   'the execution information')

    parser.add_argument('-l', '--readlen', dest="readlen", type=int, default='76', help='Length of reads in the '
                                                                                        'samfile"')
    parser.add_argument('-p', '--prefix', dest="prefix", type=str, default='', help='Output prefix string to '
                                                                                    'personalize partially the output '
                                                                                    'file.')
    parser.add_argument('-g', '--genome', dest="genome", help='Genome version an species"')
    parser.add_argument('--pcr', dest='pcr_filename', action="store", help='PCR bed file as gold_standard')
    parser.add_argument('--gff_output', dest='gff_output', action="store", help='Filename where a gff with the lsv '
                                                                                'events will be generated')

    parser.add_argument('-t', '--nthreads', dest="nthreads", type=int, default='4', help='Number of CPUs to use')
    parser.add_argument('-o', '--output', dest='output', action="store", required=True, help='casete exon list file')
    parser.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')
    parser.add_argument('--debug', type=int, default=0, help="Activate this flag for debugging purposes, activates "
                                                             "logger and jumps some processing steps.")
    parser.add_argument('--only_gather', action='store_true', dest='onlygather', default=False)

    return parser.parse_args()


#########
# MAIN  #
#########


def main(params):

    mglobals.global_conf_ini(params.conf, params)

    logger = utils.get_logger("%s/majiq.log" % mglobals.outDir, silent=params.silent, debug=params.debug)
    logger.info("")
    logger.info("Command: %s" % params)

    if not params.onlygather:
        p = Process(target=__parallel_gff3, args=(params.transcripts, params.pcr_filename, params.output))
        logger.info("... waiting gff3 parsing")
        p.start()
        p.join()
    chr_list = majiq_io.load_bin_file("%s/tmp/chromlist.pkl" % mglobals.outDir)

    if not params.onlygather:
        logger.info("Get samfiles")
        sam_list = []
        for exp_idx, exp in enumerate(mglobals.exp_list):
            samfile = "%s/%s.bam" % (mglobals.sam_dir, exp)
            if not os.path.exists(samfile):
                logger.info("Skipping %s.... not found" % samfile)
                continue
            sam_list.append(samfile)
        #majiq_io.count_mapped_reads(samfile, exp_idx)
        if len(sam_list) == 0:
            return

        if params.nthreads > 1:
            pool = Pool(processes=params.nthreads)
        logger.info("Scatter in Chromosomes")
        for chrom in chr_list:
            temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
            utils.create_if_not_exists(temp_dir)
            if params.nthreads == 1:
                majiq_builder(sam_list, chrom, pcr_validation=params.pcr_filename, gff_output=params.gff_output,
                          only_rna=params.only_rna, nondenovo=params.non_denovo, logging=logger)
            else:

                pool.apply_async(__parallel_lsv_quant, [sam_list, chrom,
                                                    params.pcr_filename,
                                                    params.gff_output,
                                                    params.only_rna,
                                                    params.non_denovo])

        if params.nthreads > 1:
            logger.info("... waiting childs")
            pool.close()
            pool.join()

    utils.gc_factor_calculation(chr_list, 10)
    utils.plot_gc_content()

    #GATHER
    logger.info("Gather outputs")
    utils.merge_and_create_majiq_file(chr_list, params.prefix)

    if not params.gff_output is None:
        logger.info("Gather lsv and generate gff")
        fp = open('%s/%s' % (mglobals.outDir, params.gff_output), 'w+')
        for chrom in chr_list:
            temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
            yfile = '%s/temp_gff.pkl' % temp_dir
            gff_list = majiq_io.load_bin_file(yfile)
            for gff in gff_list:
                fp.write("%s\n" % gff)
        fp.close()

    if not params.pcr_filename is None:
        logger.info("Gather pcr results")
        fp = open('%s/pcr_match.tab' % mglobals.outDir, 'w+')
        for chrom in chr_list:
            temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
            yfile = '%s/pcr.pkl' % temp_dir
            pcr_l = majiq_io.load_bin_file(yfile)
            for pcr in pcr_l:
                fp.write("%s\n" % pcr)
        fp.close()

    mglobals.print_numbers()
    logger.info("End of execution")

if __name__ == "__main__":
    args = _generate_parser()
    main(args)
