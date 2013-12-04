#!/usr/bin/python

#standard python libraries

import argparse,math, os
import numpy as np
import copy
import multiprocessing
from multiprocessing import Pool, Manager

try:
    import cPickle as pickle
except:
    import pickle
import grimoire.analize as analize
import grimoire.rnaseq_io as rnaseq_io
import grimoire.utils.utils as utils
import grimoire.mglobals as mglobals

def __parallel_for_splc_quant(samfiles_list, gene_list, chr, as_db):
    print "START child,", multiprocessing.current_process().name
    for idx,exp in enumerate(samfiles_list):
        print "READING ", exp
        rnaseq_io.reads_for_junc_coverage(exp, gene_list, mglobals.readLen, idx )

    analize.analize_genes(gene_list, "kk", as_db, None, 'AS')
    analize.analize_genes(gene_list, "kk", as_db, None, 'CONST')
    a,b = analize.analize_junction_reads( gene_list,chr )
    file_name = 'genelist_%s.obj'%(chr)
    utils.prepare_MAJIQ_table( a,b,file_name)
#    pickle.dump((tiss,exp,p2p), file_pi)
#    file_pi.close()
    print "END child, ", multiprocessing.current_process().name

def __parallel_for_body(SAM, all_genes, exp_idx, chr_list, read_len):
    exp = mglobals.exp_list[exp_idx]
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
    
    if os.path.exists(args.conf) :
        mglobals.global_conf_ini(args.conf)
    else:
        mglobals.global_init(args.readlen, args.dir, args.paths)

    all_genes = rnaseq_io.read_transcript_ucsc(args.transcripts)
    av_altern = rnaseq_io.read_triplets_bed("/data/ucsc/reads/test_1k/annotated_db/alt.chr1.sorted.mm10.bed", all_genes)

    n_genes = 0
    for gg in all_genes.values(): n_genes += len(gg) 

    if int(args.ncpus) >1:
        pool = Pool(processes=args.ncpus)              # start 4 worker processes
    chr_list = all_genes.keys()

    jobs = []
    for chrom in chr_list:
        sam_list = []
        for exp in mglobals.exp_list:
            SAM = "%s/%s.%s.sorted.sam"%(mglobals.sam_dir,exp,chrom)
            if not os.path.exists(SAM): continue
            sam_list.append(SAM)
        if len(sam_list) == 0: continue
        if int(args.ncpus) == 1:
            __parallel_for_splc_quant(sam_list, all_genes[chrom], chrom, av_altern)
        else:
            jobs.append(pool.apply_async(__parallel_for_splc_quant, [SAM,all_genes[chr], chr], av_altern))

    print "MASTER JOB.... waiting childs"
    genes = np.zeros(shape=(len(mglobals.exp_list)),dtype=np.dtype('object'))
    if int(args.ncpus) >1:
        pool.close()
        for idx, j in enumerate(jobs):
            gene[idx]=j.get()
            
        pool.join()

    #print "GEN MATLAB"
    #tiss = {}
    #exp = {}
    #p2p = {}
    #for chrom in chr_list:
    #    filename = 'genelist_%s.obj'%chrom
    #    if not os.path.exists(filename): continue
    #    file_pi2 = open(filename, 'rb')
    #    tiss[chrom],exp[chrom],p2p[chrom] = pickle.load(file_pi2)
    #
    #utils.merge_and_create_MAJIQ_matlab(tiss,exp,p2p)





if __name__ == "__main__":
    main()
