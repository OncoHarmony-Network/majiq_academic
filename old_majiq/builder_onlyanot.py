#!/usr/bin/python

import argparse
import sys
import traceback
from multiprocessing import current_process

import numpy as np

import majiq.grimoire.lsv as majiq_lsv
import majiq.src.io_utils
import old_majiq.grimoire.exon as majiq_exons
import old_majiq.grimoire.junction as majiq_junction
import old_majiq.src.voila_wrapper
from old_majiq.grimoire.gene import recreate_gene_tlb

# import old_majiq.src.analize as analize
import old_majiq.src.io as majiq_io
from old_majiq.src.normalize import prepare_gc_content
import majiq.src.utils as utils
import majiq.src.config as config
from majiq.grimoire.lsv import SSOURCE, STARGET, InvalidLSV

try:
    import cPickle as pickle
except Exception:
    import pickle


def lsv_detection(gene_list):
    lsv_list = {}
    const_set = {}

    for name, ind_list in config.tissue_repl.items():
        lsv_list[name] = []
        const_set[name] = set()

    for gn in gene_list:
        local_const = set(gn.get_all_junctions())
        local_lsv_jun = {}
        dummy = {}
        for name, ind_list in config.tissue_repl.items():
            dummy[name] = [[], []]
        for ex in gn.get_exon_list():
            try:
                ex.detect_lsv(gn, SSOURCE, dummy, local_lsv_jun, only_annot=True)
            except InvalidLSV as e:
                #print e.msg
                pass

            try:
                ex.detect_lsv(gn, STARGET, dummy, local_lsv_jun, only_annot=True)
            except InvalidLSV as e:
                #print e.msg
                pass

        for name, ind_list in config.tissue_repl.items():

            local_const.difference(local_lsv_jun)
            const_set[name].update(local_const)

            for ss in dummy[name][0]:
                for st in dummy[name][1]:
                    if ss.contained(st):
                        break
                else:
                    lsv_list[name].append(ss)

            for st in dummy[name][1]:
                for ss in dummy[name][0]:
                    if st.contained(ss):
                        break
                else:
                    lsv_list[name].append(st)

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
        tlogger = utils.get_logger("%s/db.old_majiq.log" % config.outDir, silent=False, debug=0)
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
    import random
    for name, ind_list in config.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):
            all_visual = []
            as_table = []
            nonas_table = []
            for chnk in range(config.num_final_chunks):
                temp_dir = "%s/tmp/chunk_%s" % (config.outDir, chnk)
                temp_filename = '%s/%s.splicegraph.pkl' % (temp_dir, config.exp_list[exp_idx])
                visual_gene_list = majiq.src.io_utils.load_bin_file(temp_filename)
                all_visual.append(visual_gene_list)
            fname = '%s/%s.splicegraph' % (config.outDir, config.exp_list[exp_idx])
            visual = np.concatenate(all_visual)
            majiq_io.dump_bin_file(visual, fname)
            del all_visual
            del visual

            for chnk in range(config.num_final_chunks):
                temp_dir = "%s/tmp/chunk_%s" % (config.outDir, chnk)
                filename = "%s/%s.old_majiq.pkl" % (temp_dir, config.exp_list[exp_idx])
                temp_table = majiq.src.io_utils.load_bin_file(filename)
                as_table.append(temp_table[0])
                nonas_table.append(temp_table[1])

            if len(as_table) == 0:
                continue

            info = dict()
            info['experiment'] = config.exp_list[exp_idx]
            info['genome'] = config.genome

            at = np.concatenate(as_table)
            for lsv in at:
                lsv.set_gc_factor(exp_idx)
            nat = np.concatenate(nonas_table)

            clist = random.sample(nat, min(5000, len(nat)))
            for jnc in clist:
                jnc.set_gc_factor(exp_idx)

            fname = '%s/%s.old_majiq' % (config.outDir, config.exp_list[exp_idx])
            majiq_io.dump_bin_file((info, at, clist), fname)


def prepare_intronic_exons(gene_list):

    for gn in gene_list:
        ex_list = gn.get_exon_list()
        for st, end in gn.get_ir_definition():

            exon1 = None
            exon2 = None
            for ex in ex_list:
                coords = ex.get_coordinates()
                if coords[1] == st - 1:
                    exon1 = ex
                elif coords[0] == end + 1:
                    exon2 = ex

            if exon1 is None or exon2 is None:
                continue

            junc1 = majiq_junction.Junction(st-1, st, exon1, None, gn, readN=0)
            junc2 = majiq_junction.Junction(end, end+1, None, exon2, gn, readN=0)

            ex = majiq_exons.Exon(st, end, gn, gn.get_strand(), annot=True, isintron=True)
            gn.add_exon(ex)

            trcpt = exon1.exonTx_list[0].get_transcript()

            txex = majiq_exons.ExonTx(st, end, trcpt[0], intron=True)
            txex.add_3prime_junc(junc1)
            txex.add_5prime_junc(junc2)

            ex.ss_3p_list.append(txex.start)
            ex.ss_5p_list.append(txex.end)
            ex.exonTx_list.append(txex)
                # txex.exon = ex
            junc1.add_acceptor(ex)
            junc2.add_donor(ex)

            junc1.add_donor(exon1)
            for ex in exon1.exonTx_list:
                st, end = ex.get_coordinates()
                if end == junc1.get_coordinates()[0]:
                    ex.add_5prime_junc(junc1)
                    break

            junc2.add_acceptor(exon2)
            for ex in exon2.exonTx_list:
                st, end = ex.get_coordinates()
                if st == junc2.get_coordinates()[1]:
                    ex.add_3prime_junc(junc2)
                    break
            break
        gn.prepare_exons()


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

    config.global_conf_ini(params.conf, params, only_db=True)

    logger = utils.get_logger("%s/old_majiq.log" % config.outDir)
    logger.info("")
    logger.info("Command: %s" % params)

    majiq_io.read_gff(params.transcripts, None, 1, logging=logger)
    chr_list = majiq.src.io_utils.load_bin_file("%s/tmp/chromlist.pkl" % config.outDir)

    for chnk in range(config.num_final_chunks):

        if not logger is None:
            logger.info("Building for chunk %s" % chnk)

        temp_dir = "%s/tmp/chunk_%s" % (config.outDir, chnk)
        temp_file = open('%s/annot_genes.pkl' % temp_dir, 'rb')
        gene_list = pickle.load(temp_file)
        prepare_intronic_exons(gene_list)
        if not logger is None:
            logger.info("[%s] Recreatin Gene TLB" % chnk)
        recreate_gene_tlb(gene_list)

        if not logger is None:
            logger.info("[%s] Detecting LSV" % chnk)
        lsv, const = lsv_detection(gene_list)

        prepare_gc_content(gene_list, temp_dir)

        old_majiq.src.voila_wrapper.generate_visualization_output(gene_list, temp_dir)
        if not logger is None:
            logger.info("[%s] Preparing output" % chnk)

        utils.prepare_lsv_table(lsv, const, temp_dir)
        majiq_lsv.extract_gff(lsv, temp_dir)
        #ANALYZE_DENOVO
        # utils.analyze_denovo_junctions(gene_list, "%s/denovo.pkl" % temp_dir)
        # utils.histogram_for_exon_analysis(gene_list, "%s/ex_lengths.pkl" % temp_dir)

    #GATHER
    logger.info("Gather outputs")
    merge_and_create_majiq_file(chr_list, config.outDir)

    logger.info("Gather lsv and generate gff")
    fp = open('%s/%s' % (config.outDir, 'lsvs.gff3'), 'w+')
    for chnk in range(config.num_final_chunks):
        temp_dir = "%s/tmp/chunk_%s" % (config.outDir, chnk)
        yfile = '%s/temp_gff.pkl' % temp_dir
        gff_list = majiq.src.io_utils.load_bin_file(yfile)
        for gff in gff_list:
            fp.write("%s\n" % gff)
    fp.close()

    config.print_numbers()
    logger.info("End of execution")

if __name__ == "__main__":
    args = _generate_parser()
    main(args)
