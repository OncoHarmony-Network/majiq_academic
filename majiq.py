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

def __parallel_for_splc_quant(samfiles_list, gene_list, chr, as_db):
    print "START child,", current_process().name
    for idx,exp in enumerate(samfiles_list):
        print "READING ", idx, exp
        rnaseq_io.read_sam_or_bam(exp, gene_list, mglobals.readLen, chr, idx )
#        rnaseq_io.reads_for_junc_coverage(exp, gene_list, mglobals.readLen, idx )
    analize.annotated_AS_events(gene_list, 'AS')
    a,b = analize.rnaseq_AS_events( gene_list, chr )
    file_name = '%s.obj'%(chr)
    utils.prepare_MAJIQ_table( a,b,file_name)

    print "END child, ", current_process().name

def __parallel_for_body(SAM, all_genes, exp_idx, chr_list, read_len):
    exp = mglobals.exp_list[exp_idx]
    print "START child,",current_process().name
    rnaseq_io.reads_for_junc_coverage(SAM,all_genes,read_len,exp_idx)
    print "END child, ", current_process().name

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
    parser.add_argument('-t','--ncpus', dest="ncpus", type=int,default='4', help='Number of CPUs to use')
    parser.add_argument('-o','--output', dest='output',action= "store", help='casete exon list file')
    return parser.parse_args()


#########
# MAIN  #
#########

def main() :

    args = _generate_parser()
    print args
    if os.path.exists(args.conf) :
        mglobals.global_conf_ini(args.conf)
    else:
        print "NO CONF"
        mglobals.global_init(args.readlen, args.dir, args.paths)

    all_genes = rnaseq_io.read_transcript_ucsc(args.transcripts)
#    av_altern = rnaseq_io.read_triplets_bed("/data/ucsc/reads/test_1k/annotated_db/alt.chr1.sorted.mm10.bed", all_genes)

    n_genes = 0
    temp = []
    for gl in all_genes.values(): 
        n_genes += len(gl) 
#        for genes_l in gl.values():
#            for gg in genes_l:
#                temp += gg.get_exon_list()

    if int(args.ncpus) >1:
        pool = Pool(processes=args.ncpus)              # start 4 worker processes
    chr_list = all_genes.keys()

    jobs = []
    sam_list = []
    for exp_idx, exp in enumerate(mglobals.exp_list):
#        SAM = "%s/%s.%s.sorted.sam"%(mglobals.sam_dir,exp,chrom)
        SAM = "%s/%s.sorted.bam"%(mglobals.sam_dir,exp)
        print SAM
        if not os.path.exists(SAM): continue
        sam_list.append(SAM)
        rnaseq_io.count_mapped_reads(SAM,exp_idx)
    if len(sam_list) == 0: return



    for chrom in chr_list:
        if int(args.ncpus) == 1:
            __parallel_for_splc_quant(sam_list, all_genes[chrom], chrom, None)
        else:
            jobs.append(pool.apply_async(__parallel_for_splc_quant, [sam_list, all_genes[chr], chr, av_altern]))



    
    print "MASTER JOB.... waiting childs"
    genes = np.zeros(shape=(len(mglobals.exp_list)),dtype=np.dtype('object'))
    if int(args.ncpus) >1:
        pool.close()
        for idx, j in enumerate(jobs):
            gene[idx]=j.get()
            
        pool.join()

    for gl in all_genes.values(): 
        for genes_l in gl.values():
            for gg in genes_l:
                temp += gg.get_exon_list()


    print "number of gcs", len(temp)
    utils.gc_factor_calculation(temp, 10)
    utils.plot_gc_content()
    #print "GEN MATLAB"
   # tiss = {}
   # exp = {}
   # p2p = {}

    utils.merge_and_create_MAJIQ( chr_list, 'tojuan.majiq')

if __name__ == "__main__":
    main()
