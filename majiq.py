#!/usr/bin/python

#standard python libraries

import argparse,math, os
import numpy as np
import copy
from multiprocessing import Pool, Manager, current_process
import grimoire.analize as analize
import grimoire.rnaseq_io as rnaseq_io
import grimoire.utils.utils as utils
import grimoire.mglobals as mglobals

try:
    import cPickle as pickle
except:
    import pickle

def __parallel_for_splc_quant(samfiles_list, gene_list, chr, as_db, pcr_validation = {}):
    print "START child,", current_process().name
    for idx,exp in enumerate(samfiles_list):
        print "READING ", idx, exp
        rnaseq_io.read_sam_or_bam(exp, gene_list, mglobals.readLen, chr, idx )
#        rnaseq_io.reads_for_junc_coverage(exp, gene_list, mglobals.readLen, idx )
    analize.annotated_AS_events(gene_list, 'AS')
    a,b = analize.rnaseq_AS_events( gene_list, chr )
    if pcr_validation is not None:
        utils.get_validated_pcr_events(pcr_validation,a)
    file_name = '%s.obj'%(chr)
    utils.prepare_MAJIQ_table( a,b,file_name)
    print "END child, ", current_process().name

def __parallel_lsv_quant(samfiles_list, gene_list, chr, as_db):
    #print "START child,", current_process().name
    for idx,exp in enumerate(samfiles_list):
        print "READING ", idx, exp
        rnaseq_io.read_sam_or_bam(exp, gene_list, mglobals.readLen, chr, idx )
    lsv, const = analize.LSV_detection( gene_list, chr )
    file_name = '%s.obj'%(chr)
    utils.prepare_LSV_table( lsv,const ,file_name)
    #print "END child, ", current_process().name

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

def main( args ) :

    if os.path.exists(args.conf) :
        mglobals.global_conf_ini(args.conf)
    else:
        mglobals.global_init(args.readlen, args.dir, args.paths)

    all_genes = rnaseq_io.read_transcript_ucsc(args.transcripts)
    temp = []
    for gl in all_genes.values(): 
        for genes_l in gl.values():
            for gg in genes_l:
                temp += gg.get_exon_list()

    if int(args.ncpus) >1:
        pool = Pool(processes=args.ncpus)              # start 4 worker processes
    chr_list = all_genes.keys()

    jobs = []
    sam_list = []
    if args.pcr_filename is not None:
        rnaseq_io.read_bed_pcr( args.pcr_filename , all_genes)
    for exp_idx, exp in enumerate(mglobals.exp_list):
        SAM = "%s/%s.sorted.bam"%(mglobals.sam_dir,exp)
        if not os.path.exists(SAM): 
            print "Skipping %s.... notfound"%SAM
            continue
        sam_list.append(SAM)
        rnaseq_io.count_mapped_reads(SAM,exp_idx)
    if len(sam_list) == 0: return

    if args.lsv: exec_pipe = __parallel_lsv_quant
    else: exec_pipe = __parallel_for_splc_quant

    for chrom in chr_list:
        if int(args.ncpus) == 1:
            exec_pipe(sam_list, all_genes[chrom], chrom, None)
        else:
            jobs.append(pool.apply_async( exec_pipe, [sam_list, all_genes[chrom], chrom, None]))

    print "MASTER JOB.... waiting childs"
    genes = np.zeros(shape=(len(mglobals.exp_list)),dtype=np.dtype('object'))
    if int(args.ncpus) >1:
        pool.close()
        for idx, j in enumerate(jobs):
            p=j.get()
            print p
            
        pool.join()

    #for gl in all_genes.values(): 
    #    for genes_l in gl.values():
    #        for gg in genes_l:
    #            temp += gg.get_exon_list()

    utils.generate_visualization_output(all_genes)
    print "number of gcs", len(temp)
    utils.gc_factor_calculation(temp, 10)
    utils.plot_gc_content()
    utils.merge_and_create_MAJIQ( chr_list, 'tojuan.majiq')
    mglobals.print_numbers()


if __name__ == "__main__":
    args = _generate_parser()
    print args
    main( args )
