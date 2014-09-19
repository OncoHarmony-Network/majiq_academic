#!/usr/bin/python

#standard python libraries

import argparse
import os
import sys
from multiprocessing import Pool, current_process
import grimoire.analize as analize
import grimoire.rnaseq_io as rnaseq_io
import grimoire.utils.utils as utils
import grimoire.mglobals as mglobals
import grimoire.lsv as majiq_lsv

try:
    import cPickle as pickle
except Exception:
    import pickle


def majiq_builder(samfiles_list, chrom, pcr_validation=False, logging=None):

    logging.info("Building for chromosome %s" % chrom)
    temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
    temp_file = open('%s/annot_genes.pkl' % temp_dir, 'rb')
    gene_list = pickle.load(temp_file)

    utils.create_if_not_exists(temp_dir)
    logging.info("[%s] Reading BAM files" % chrom)
    rnaseq_io.read_sam_or_bam(samfiles_list, gene_list, mglobals.readLen, chrom, logging=logging)
    logging.info("[%s] Detecting LSV" % chrom)
    lsv, const = analize.lsv_detection(gene_list, chrom, logging=logging)
    file_name = '%s.obj' % chrom
    if pcr_validation:
        utils.get_validated_pcr_lsv(lsv, temp_dir)

    majiq_lsv.extract_gff(lsv, temp_dir)
    #utils.generate_visualization_output(gene_list)
    logging.info("[%s] Preparing output" % chrom)
    utils.prepare_LSV_table(lsv, const, file_name)


def __parallel_lsv_quant(samfiles_list, gene_list, chrom, pcr_validation=False):

    try:
        print "START child,", current_process().name
        tlogger = utils.get_logger("%s/majiq.log" % mglobals.outDir, silent=args.silent, debug=args.debug)
        majiq_builder(samfiles_list, gene_list, chrom, pcr_validation, tlogger)
        print "END child, ", current_process().name
    except Exception as e:
        print "Line %s:" % sys.exc_traceback.tb_lineno, e
        sys.stdout.flush()
        raise()


def _new_subparser():
    return argparse.ArgumentParser(add_help=False)


def _generate_parser():
    dir_or_paths = _new_subparser()
    mutexc = dir_or_paths.add_mutually_exclusive_group(required=True)
    mutexc.add_argument('-dir', action="store", help='Provide a directory with all the files')
    mutexc.add_argument('-paths', default=None, nargs='+', help='Provide the list of files to analyze')
    mutexc.add_argument('-conf', default=None, help='Provide study configuration file with all '
                                                    'the execution information')

    parser = argparse.ArgumentParser(parents=[dir_or_paths])
    parser.add_argument('transcripts', action="store", help='read file in SAM format')
    parser.add_argument('-l', '--readlen', dest="readlen", type=int, default='76', help='Length of reads in the '
                                                                                        'samfile"')
    parser.add_argument('-g', '--genome', dest="genome", help='Genome version an species"')
    parser.add_argument('-pcr', dest='pcr_filename', action="store", help='PCR bed file as gold_standard')
    parser.add_argument('-t', '--ncpus', dest="ncpus", type=int, default='4', help='Number of CPUs to use')
    parser.add_argument('-o', '--output', dest='output', action="store", help='casete exon list file')
    parser.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')
    parser.add_argument('--debug', type=int, default=0, help="Activate this flag for debugging purposes, activates "
                                                             "logger and jumps some processing steps.")
    return parser.parse_args()

#########
# MAIN  #
#########


def main(params):

    if os.path.exists(params.conf):
        mglobals.global_conf_ini(params.conf)
    else:
        mglobals.global_init(params.readlen, params.dir, params.paths)

    mglobals.gene_tlb = {}

    chr_list = rnaseq_io.read_gff(params.transcripts)
    #chr_list = all_genes.keys()
    # temp = []
    # for gl in all_genes.values():
    #     for genes_l in gl.values():
    #         for gg in genes_l:
    #             temp += gg.get_exon_list()

    # if params.pcr_filename is not None:
    #     rnaseq_io.read_bed_pcr(params.pcr_filename, all_genes)

    logger = utils.get_logger("%s/majiq.log" % mglobals.outDir, silent=args.silent, debug=args.debug)
    logger.info("")
    logger.info("Command: %s" % params)
    sam_list = []
    for exp_idx, exp in enumerate(mglobals.exp_list):
        samfile = "%s/%s.sorted.bam" % (mglobals.sam_dir, exp)
        if not os.path.exists(samfile):
            logger.info("Skipping %s.... not found" % samfile)
            continue
        sam_list.append(samfile)
        rnaseq_io.count_mapped_reads(samfile, exp_idx)
    if len(sam_list) == 0:
        return

    if int(params.ncpus) > 1:
        pool = Pool(processes=params.ncpus)

    for chrom in chr_list:
        if int(params.ncpus) == 1:
            majiq_builder(sam_list, chrom, pcr_validation=params.pcr_filename, logging=logger)
        else:
            pool.apply_async(majiq_builder, [sam_list, chrom, params.pcr_filename])

    if int(params.ncpus) > 1:
        logger.info("... waiting childs")
        pool.close()
        pool.join()

    # utils.gc_factor_calculation(temp, 10)
    #utils.plot_gc_content()

    #GATHER
    utils.merge_and_create_majiq_file(chr_list, 'tojuan.majiq')

    fp = open('%s/lsv_miso.gtf' % mglobals.outDir, 'w+')
    fp2 = open('%s/pcr_match.tab' % mglobals.outDir, 'w+')
    for chrom in chr_list:
        temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
        yfile = '%s/temp_gff.pkl' % temp_dir
        if not os.path.exists(yfile):
            continue
        gtf_list = pickle.load(open(yfile, 'rb'))
        for gtf in gtf_list:
            fp.write("%s\n" % gtf)
        yfile = '%s/pcr.pkl' % temp_dir
        if not os.path.exists(yfile):
            continue
        pcr_l = pickle.load(open(yfile, 'rb'))
        for pcr in pcr_l:
            fp2.write("%s\n" % pcr)
    fp2.close()
    fp.close()

    mglobals.print_numbers()


if __name__ == "__main__":
    args = _generate_parser()
    main(args)
