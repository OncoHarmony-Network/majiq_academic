#!/usr/bin/python

#standard python libraries

import argparse, os
import numpy as np

import sys
from multiprocessing import Pool, current_process
import grimoire.analize as analize
import grimoire.rnaseq_io as rnaseq_io
import grimoire.utils.utils as utils
import grimoire.mglobals as mglobals
import grimoire.lsv  as majiq_lsv

try:
    import cPickle as pickle
except:
    import pickle


def majiq_builder(samfiles_list, gene_list, chrom, as_db, pcr_validation=False):

    temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
    utils.create_if_not_exists(temp_dir)
    rnaseq_io.read_sam_or_bam(samfiles_list, gene_list, mglobals.readLen, chrom)
    lsv, const = analize.LSV_detection(gene_list, chrom)
    file_name = '%s.obj' % chrom
    if pcr_validation:
        utils.get_validated_pcr_lsv(lsv, temp_dir)

    majiq_lsv.extract_gff(lsv, temp_dir)
    utils.prepare_LSV_table(lsv, const, file_name)



def __parallel_lsv_quant(samfiles_list, gene_list, chrom, as_db, pcr_validation=False):

    try:
        print "START child,", current_process().name
        majiq_builder(samfiles_list, gene_list, chrom, as_db, pcr_validation)
        print "END child, ", current_process().name
    except Exception as e:
        print "Line %s:"%sys.exc_traceback.tb_lineno, e1
        sys.stdout.flush()
        raise()

def _new_subparser():
    return argparse.ArgumentParser(add_help=False)

def _generate_parser():
    dir_or_paths = _new_subparser()
    mutexc = dir_or_paths.add_mutually_exclusive_group(required=True)
    mutexc.add_argument('-dir', action= "store", help='Provide a directory with all the files')
    mutexc.add_argument('-paths', default=None, nargs='+', help='Provide the list of files to analyze')
    mutexc.add_argument('-conf',default=None, help='Provide study configuration file with all the execution information')

    parser = argparse.ArgumentParser(parents=[dir_or_paths])
    parser.add_argument('transcripts', action= "store", help='read file in SAM format')
    parser.add_argument('-l','--readlen', dest="readlen", type=int,default='76', help='Length of reads in the samfile"')
    parser.add_argument('-g','--genome', dest="genome", help='Genome version an species"')
    parser.add_argument('-pcr', dest='pcr_filename',action= "store", help='PCR bed file as gold_standard')
    parser.add_argument('-lsv', dest="lsv", action="store_true", default=False, help='Using lsv analysis')
    parser.add_argument('-t','--ncpus', dest="ncpus", type=int,default='4', help='Number of CPUs to use')
    parser.add_argument('-o','--output', dest='output',action= "store", help='casete exon list file')
    return parser.parse_args()

#########
# MAIN  #
#########


def main(params):

    if os.path.exists(params.conf):
        mglobals.global_conf_ini(params.conf)
    else:
        mglobals.global_init(params.readlen, params.dir, params.paths)

    all_genes = rnaseq_io.read_gff(params.transcripts)
    temp = []
    for gl in all_genes.values(): 
        for genes_l in gl.values():
            for gg in genes_l:
                temp += gg.get_exon_list()

    chr_list = all_genes.keys()

    if params.pcr_filename is not None:
        rnaseq_io.read_bed_pcr(params.pcr_filename, all_genes)

    sam_list = []
    for exp_idx, exp in enumerate(mglobals.exp_list):
        samfile = "%s/%s.sorted.bam" % (mglobals.sam_dir, exp)
        if not os.path.exists(samfile):
            print "Skipping %s.... not found" % samfile
            continue
        sam_list.append(samfile)
        rnaseq_io.count_mapped_reads(samfile, exp_idx)
    if len(sam_list) == 0:
        return

    if int(params.ncpus) > 1:
        pool = Pool(processes=params.ncpus)

    for chrom in chr_list:
        if int(params.ncpus) == 1:
            majiq_builder(sam_list, all_genes[chrom], chrom, None, pcr_validation=params.pcr_filename)
        else:
            pool.apply_async(majiq_builder, [sam_list, all_genes[chrom], chrom, None, params.pcr_filename])

    print "MASTER JOB.... waiting childs"
    if int(params.ncpus) > 1:
        pool.close()
        pool.join()

    utils.generate_visualization_output(all_genes)
    print "number of gcs", len(temp)
    utils.gc_factor_calculation(temp, 10)
    #utils.plot_gc_content()
    utils.merge_and_create_MAJIQ( chr_list, 'tojuan.majiq')

    fp = open('%s/lsv_miso.gtf' % (mglobals.outDir), 'w+')
    fp2 = open('%s/pcr_match.tab' % (mglobals.outDir), 'w+')
    for chrom in chr_list:
        temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
        yfile = '%s/temp_gff.pkl' % (temp_dir)
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
    print args
    main(args)
