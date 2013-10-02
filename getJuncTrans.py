#!/usr/bin/python

#from __future__ import division
import argparse,math, os
import numpy as np
import analize
import globals 
import rnaseq_io
from utils import utils
import copy
#from matplotlib import pyplot as plt
import multiprocessing
from multiprocessing import Pool

try:
    import cPickle as pickle
except:
    import pickle


def __parallel_for_splc_quant(sam_dir, gene_list,n_genes, chr, order,as_db):
#    reference_genes = copy.deepcopy(all_genes)
    print "START child,", multiprocessing.current_process().name
    for idx,exp in enumerate(globals.exp_list):
        SAM = sam_dir+"/"+exp+".sorted.sam"
        rnaseq_io.reads_for_junc_coverage(SAM, gene_list, globals.readLen, idx )
    print "END child, ", multiprocessing.current_process().name
    TAS = analize.analize_genes(gene_list, "kk", as_db, None, 'AS')
    CONST = analize.analize_genes(gene_list, "kk", as_db, None, 'CONST')
    a,b = analize.analize_junction_reads( gene_list,chr )
    return (a,b)

def __parallel_for_body(SAM, all_genes,n_genes, exp_idx,chr_list, order, read_len):
#    reference_genes = copy.deepcopy(all_genes)
    exp = globals.exp_list[exp_idx]
    print "START child,", multiprocessing.current_process().name
#    gen = np.zeros( shape=(n_genes), dtype=np.dtype('object'))
    rnaseq_io.reads_for_junc_coverage(SAM,all_genes,read_len,exp_idx)


#    print "Creating matlab file"
#    utils.create_junction_matlab_file(all_genes,'./test_'+exp+'.mat', chr_list,order)
    print "END child, ", multiprocessing.current_process().name
    return
#    return gen


# MAIN

if __name__ == "__main__":
    #parser for the basic flags that all subcommands share
    basic_flags = argparse.ArgumentParser(add_help=False)
    basic_flags.add_argument('-l','--readlen', dest="readlen", type=int,default='76', help='Length of reads in the samfile"')
    basic_flags.add_argument('-g','--genome', dest="genome", help='Genome version an species"')
    basic_flags.add_argument('-t','--ncpus', dest="ncpus", type=int,default='4', help='Number of CPUs to use')
    #main parser    
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='operation')
    #subparser "input", to start from the SAM files
    parser_input = subparsers.add_parser('input', help='Input file and create binary file', parents=[basic_flags])
    parser_input.add_argument('transcripts', action= "store", help='read file in SAM format')
    parser_input.add_argument('junc_dir', action= "store", help='reads directory in SAM format')
    parser_input.add_argument('-r', action='store_true', default=False)
    parser_input.add_argument('-o','--output', dest='output',action= "store", help='casete exon list file')
    # subparser "load", to start from processed files
    parser_load = subparsers.add_parser('load', help='Load binary file', parents=[basic_flags])
    parser_load.add_argument('bin', action= "store", help='read binary file')
    parser_load.add_argument('junc', action= "store", help='read read binary juntion file')
    parser_load.add_argument('-o','--output', dest='output',action= "store", help='casete exon list file')
    # subparser "test"
    parser_input = subparsers.add_parser('test', help='Execute test file', parents=[basic_flags])
    parser_input.add_argument('test_file', action= "store", help='read file in SAM format')
    
    args = parser.parse_args()
    print args

    if args.operation == 'test':
        read_test_file(args.transcripts)
       
    elif args.operation == 'input':
        globals.global_init(args.readlen, args.junc_dir)
        print globals.num_experiments,globals.readLen
        all_genes = rnaseq_io.read_transcript_ucsc(args.transcripts,refSeq=args.r)
        av_altern = rnaseq_io.read_triplets_bed("/data/ucsc/reads/test_1k/annotated_db/alt.chr1.sorted.mm10.bed",all_genes)
        print "AV_ALTERN:",len(av_altern)
        n_genes = 0
        for chr,gg in all_genes.items():
            n_genes += len(gg) 
        print n_genes
        pool = Pool(processes=2)              # start 4 worker processes
        jobs = []
        order= {}
        chr_list = all_genes.keys()
#        for chr in chr_list:
#            if not chr in order:
#                order[chr]=[]
#            order[chr] = all_genes[chr]

#        print TAS
        cand = {}
        non_cand = {}
        for chr in chr_list:
#            all_experiments[idx] = RNAexperiment(exp,1,args.genome,all_genes)
#            gene_list = all_experiments[idx].get_gene_list()
            if int(args.ncpus) == 1:
                cand[chr],non_cand[chr] = __parallel_for_splc_quant(args.junc_dir,all_genes[chr],n_genes,chr,order,av_altern)
#                __parallel_for_body(SAM,all_genes,n_genes,idx,chr_list,order,args.readlen)
            else:
                jobs.append(pool.apply_async(__parallel_for_splc_quant, [SAM,all_genes[chr],n_genes,chr,order],av_alter ))
#                jobs.append(pool.apply_async(__parallel_for_body,[SAM, all_genes, n_genes, idx, chr_list, order, args.readlen] ))

        print "MASTER JOB.... waiting childs"
        genes = np.zeros(shape=(len(globals.exp_list)),dtype=np.dtype('object'))
        if int(args.ncpus) >1:
            pool.close()
            for idx, j in enumerate(jobs):
                gene[idx]=j.get()
            pool.join()

        print "GEN MATLAB"

        utils.prepare_MAJIQ_matlab_table(cand,non_cand)


#            print idx,command
#            if output: print output
#        print "read alt"
#        print "read const"
#        av_const = read_triplets_bed("/data/ucsc/const.sorted.mm10.bed", all_genes)
#        fpc = open("./wei.coord","w+")
#        for ii in av_const:
#            ii[1].print_triplet_coord(fpc)
#        fpc.close()
#        
#        print "Analize"
##        TCONS = analize.analize_genes(all_genes, args.output, av_const,av_altern,'Const')
##        TAS = analize.analize_genes(all_genes, args.output, av_altern,av_const,'AS')

#        get_junctions_STAR(args.junc,all_genes)
    

     #prepare_const_set( av_const, TCONS, all_rest[1] , user_lab='CONST')
#    prepare_const_set( av_altern, TAS, all_rest[0], user_lab = 'AS')
#    I = analize_junction_reads(all_genes, av_altern,av_const, (False,True,True) )
#    II = analize_junction_reads(all_genes, av_altern,av_const, (True,False,True) )
#    III = analize_junction_reads(all_genes, av_altern,av_const, (True,True,False) )
#    venn_gen(I,II,III)
