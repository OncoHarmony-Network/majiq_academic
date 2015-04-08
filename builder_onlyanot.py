#!/usr/bin/python

import argparse
import os
import sys
import traceback
from multiprocessing import Pool, current_process, Process
import grimoire.analize as analize
import grimoire.rnaseq_io as majiq_io
import grimoire.utils.utils as utils
import grimoire.mglobals as mglobals
import numpy as np
import grimoire.lsv as majiq_lsv

try:
    import cPickle as pickle
except Exception:
    import pickle


def lsv_detection(gene_list, chrom, only_real_data=False, logging=None):

    num_ss_var = [[0]*20, [0]*20, 0]

    const_set = [set() for xx in range(mglobals.num_experiments)]
    lsv_list = [[] for xx in range(mglobals.num_experiments)]
    jun = {}
    for xx in mglobals.tissue_repl.keys():
        jun[xx] = set()

    for strand, glist in gene_list.items():
        for gn in glist:

            gn.check_exons()
            mat, exon_list, tlb, var_ss = gn.get_rnaseq_mat(const_set, use_annot=not only_real_data)
            vip = []
            for idx, ex in enumerate(exon_list):
                sc = ex.get_pcr_score()
                if sc is None:
                    continue
                vip.append(idx)

            for ss in range(2):
                for ssnum in range(20):
                    num_ss_var[ss][ssnum] += var_ss[ss][ssnum]
            num_ss_var[2] += var_ss[2]
            #num_ss_var [1]+= var_ss[1]

#            print "---------------- %s --------------"%gn.get_id()
#             utils.print_junc_matrices(mat, tlb=tlb, fp=True)
            SS, ST = lsv_matrix_detection(mat, tlb, (False, False, False), vip)
            dummy = {}
            for name, ind_list in mglobals.tissue_repl.items():
                dummy[name] = [[], []]

            for lsv_index, lsv_lst in enumerate((SS, ST)):
                lsv_type = (analize.SSOURCE, analize.STARGET)[lsv_index]
                sstype = ['5prime', '3prime'][lsv_index]
#                print lsv_lst

                for idx in lsv_lst:
                    jlist = exon_list[idx].get_junctions(sstype)
                    jlist = [x for x in jlist if x is not None]
                    if len(jlist) == 0:
                        continue

                    lsv_in = gn.new_lsv_definition(exon_list[idx], jlist, lsv_type)
                    if lsv_in is None:
                        continue

                    for name, ind_list in mglobals.tissue_repl.items():
                        for jj in jlist:
                            jun[name].add(jj)

                        dummy[name][lsv_index].append(lsv_in)

            for name, ind_list in mglobals.tissue_repl.items():
                for ss in dummy[name][0]:
                    for st in dummy[name][1]:
                        if ss.contained(st):
                            break
                    else:
                        for exp_idx in ind_list:
                            lsv_list[exp_idx].append(ss)

                for st in dummy[name][1]:
                    for ss in dummy[name][0]:
                        if st.contained(ss):
                            break
                    else:
                        for exp_idx in ind_list:
                            lsv_list[exp_idx].append(st)

    for name, ind_list in mglobals.tissue_repl.items():
        for exp_idx in ind_list:
            const_set[exp_idx].difference(jun[name])

    return lsv_list, const_set


def lsv_matrix_detection(mat, exon_to_ss, b_list, vip_set=[]):
    """
       Rules for const are:
        1. All the junction from A should go to C1 or C2
        2. All the junction from C1 should go to A
        3. All the junction to C2 should come from A
        4. Number of reads from C1-A should be equivalent to number of reads from A-C2
    """
    lsv_list = [[], []]

    #change bucle for iterate by exons
    for ii in range(0, len(exon_to_ss) - 1):
        lsv = exon_to_ss[ii]
        #Single Source detection
        ss = mat[lsv[1][0]:lsv[1][-1]+1, :]
        ss_valid = True

        if ss_valid and np.count_nonzero(ss) > 1:
            lsv_list[0].append(ii)

    for ii in range(1, len(exon_to_ss)):
        lsv = exon_to_ss[ii]
        #Single Targe detection
        st = mat[:, lsv[0][0]:lsv[0][-1]+1]
        st_valid = True

        if st_valid and np.count_nonzero(st) > 1:
            lsv_list[1].append(ii)

    return lsv_list




def __parallel_gff3(transcripts):

    try:
        print "START child,", current_process().name
        tlogger = utils.get_logger("%s/db.majiq.log" % mglobals.outDir, silent=False, debug=0)
        majiq_io.read_gff(transcripts, None, logging=tlogger)
        print "END child, ", current_process().name
    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


def merge_and_create_majiq_file(chr_list, pref_file):

    """

    :param chr_list:
    :param pref_file:
    """
    if pref_file != '':
        pref_file = '%s.' % pref_file
    import random
    for name, ind_list in mglobals.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):
            all_visual = []
            as_table = []
            nonas_table = []
            for chrom in chr_list:
                temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
                temp_filename = '%s/%s.splicegraph.pkl' % (temp_dir, mglobals.exp_list[exp_idx])
                visual_gene_list = majiq_io.load_bin_file(temp_filename)
                all_visual.append(visual_gene_list)
            fname = '%s/%s%s.splicegraph' % (mglobals.outDir, pref_file, mglobals.exp_list[exp_idx])
            visual = np.concatenate(all_visual)
            majiq_io.dump_bin_file(visual, fname)
            del all_visual
            del visual

            for chrom in chr_list:
                temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
                filename = "%s/%s.majiq.pkl" % (temp_dir, mglobals.exp_list[exp_idx])
                temp_table = majiq_io.load_bin_file(filename)
                as_table.append(temp_table[0])
                nonas_table.append(temp_table[1])

            if len(as_table) == 0:
                continue

            info = dict()
            info['experiment'] = mglobals.exp_list[exp_idx]

            info['genome'] = mglobals.genome
            info['num_reads'] = mglobals.num_mapped_reads[exp_idx]
            at = np.concatenate(as_table)
            nat = np.concatenate(nonas_table)

            clist = random.sample(nat, min(5000, len(nat)))

            fname = '%s/%s%s.majiq' % (mglobals.outDir, pref_file, mglobals.exp_list[exp_idx])
            majiq_io.dump_bin_file((info, at, clist), fname)



def _new_subparser():
    return argparse.ArgumentParser(add_help=False)


def _generate_parser():
    parser = argparse.ArgumentParser(description="MAJIQ is a suite of tools for the analysis of Alternative "
                                                 "Splicing Events and Alternative Splicing Quantification.")
    #common flags (first ones are required)
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    parser.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    parser.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')

    parser.add_argument('transcripts', action="store", help='read file in SAM format')
    parser.add_argument('-conf', default=None, help='Provide study configuration file with all '
                                                    'the execution information')
    parser.add_argument('--minreads', default=2, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'the LSV "exist in the data". '
                             '[Default: %(default)s]')
    parser.add_argument('--minpos', default=2, type=int, help='Minimum number of start positions with at least 1 '
                                                                   'read in a LSV to consider that the LSV "exist in '
                                                                   'the data"')

    return parser.parse_args()


#########
# MAIN  #
#########


def main(params):

    mglobals.global_conf_ini(params.conf, params, only_db=True)

    logger = utils.get_logger("%s/majiq.log" % mglobals.outDir)
    logger.info("")
    logger.info("Command: %s" % params)

    majiq_io.read_gff(params.transcripts, None, logging=logger)
    chr_list = majiq_io.load_bin_file("%s/tmp/chromlist.pkl" % mglobals.outDir)

    for chrom in chr_list:

        if not logger is None:
            logger.info("Building for chromosome %s" % chrom)

        temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
        temp_file = open('%s/annot_genes.pkl' % temp_dir, 'rb')
        gene_list = pickle.load(temp_file)


        if not logger is None:
            logger.info("[%s] Recreatin Gene TLB" % chrom)
        utils.recreate_gene_tlb(gene_list)

        if not logger is None:
            logger.info("[%s] Detecting LSV" % chrom)
        lsv, const = lsv_detection(gene_list, chrom, logging=logger)

        utils.prepare_gc_content(gene_list, temp_dir)

        utils.generate_visualization_output(gene_list, temp_dir)
        if not logger is None:
            logger.info("[%s] Preparing output" % chrom)

        utils.prepare_lsv_table(lsv, const, temp_dir)

        #ANALYZE_DENOVO
        utils.analyze_denovo_junctions(gene_list, "%s/denovo.pkl" % temp_dir)
        utils.histogram_for_exon_analysis(gene_list, "%s/ex_lengths.pkl" % temp_dir)

    #GATHER
    logger.info("Gather outputs")
    merge_and_create_majiq_file(chr_list, './output')

    mglobals.print_numbers()
    logger.info("End of execution")

if __name__ == "__main__":
    args = _generate_parser()
    main(args)
