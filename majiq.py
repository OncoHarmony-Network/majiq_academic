#!/usr/bin/python

#standard python libraries
import argparse
import math
import os
import numpy as np
import copy
import multiprocessing
from multiprocessing import Pool
try:
    import cPickle as pickle
except:
    import pickle
    
#majiq libraries
import my_globals
import analize
import rnaseq_io
from utils import utils


def __parallel_for_splc_quant(sam_dir, gene_list, n_genes, chr, order,as_db):
    print "START child,", multiprocessing.current_process().name
    for path_count, path in enumerate(my_globals.exp_list):
        rnaseq_io.reads_for_junc_coverage(path, gene_list, my_globals.readLen, path_count)
        
    print "END child, ", multiprocessing.current_process().name
    TAS = analize.analize_genes(gene_list, "kk", as_db, None, 'AS')
    CONST = analize.analize_genes(gene_list, "kk", as_db, None, 'CONST')
    a,b = analize.analize_junction_reads( gene_list,chr )
    return (a,b)

def __parallel_for_body(SAM, all_genes, n_genes, exp_idx, chr_list, order, read_len):
    exp = my_globals.exp_list[exp_idx]
    print "START child,", multiprocessing.current_process().name
    rnaseq_io.reads_for_junc_coverage(SAM,all_genes,read_len,exp_idx)
    print "END child, ", multiprocessing.current_process().name

def _new_subparser():
    return argparse.ArgumentParser(add_help=False)

def _generate_parser():
    dir_or_paths = _new_subparser()
    mutexc = dir_or_paths.add_mutually_exclusive_group(required=True)
    mutexc.add_argument('-dir', action= "store", help='Provide a directory with all the files')
    mutexc.add_argument('-paths', default=None, nargs='+', help='Provide the list of files to analyze')
    
    parser = argparse.ArgumentParser(parents=[dir_or_paths])
    parser.add_argument('transcripts', action= "store", help='read file in SAM format')
    parser.add_argument('-r', action='store_true', default=False)
    parser.add_argument('-o','--output', dest='output',action= "store", help='casete exon list file')
    parser.add_argument('-l','--readlen', dest="readlen", type=int,default='76', help='Length of reads in the samfile"')
    parser.add_argument('-g','--genome', dest="genome", help='Genome version an species"')
    parser.add_argument('-t','--ncpus', dest="ncpus", type=int,default='4', help='Number of CPUs to use')
    return parser.parse_args()
                                                                                
#########
# MAIN  #
#########
if __name__ == "__main__":
    args = _generate_parser()
    my_globals.global_init(args.readlen, args.dir, args.paths)
    all_genes = rnaseq_io.read_transcript_ucsc(args.transcripts, refSeq=args.r)
    av_altern = rnaseq_io.read_triplets_bed("/data/ucsc/reads/test_1k/annotated_db/alt.chr1.sorted.mm10.bed", all_genes)
    n_genes = 0
    for gg in all_genes.values():
        n_genes += len(gg) 

    print "NUM GENES:", n_genes
    pool = Pool(processes=2)              # start 4 worker processes
    jobs = []
    order= {}
    chr_list = all_genes.keys()
    cand = {}
    non_cand = {}
    for chrom in chr_list:
        if int(args.ncpus) == 1:
            cand[chrom],non_cand[chrom] = __parallel_for_splc_quant(args.dir, all_genes[chrom], n_genes, chrom, order, av_altern)
        else:
            jobs.append(pool.apply_async(__parallel_for_splc_quant, [SAM,all_genes[chr],n_genes, chr, order], av_altern))

    print "MASTER JOB.... waiting childs"
    genes = np.zeros(shape=(len(my_globals.exp_list)),dtype=np.dtype('object'))
    if int(args.ncpus) >1:
        pool.close()
        for idx, j in enumerate(jobs):
            gene[idx]=j.get()
            
        pool.join()

    print "GEN MATLAB"
    utils.prepare_MAJIQ_matlab_table(cand,non_cand)

