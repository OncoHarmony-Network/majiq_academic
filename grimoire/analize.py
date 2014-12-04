import numpy as np
from lsv import SSOURCE, STARGET
import mglobals


def reliable_in_data(junc, exp_idx):
    
    min_read_x_exp = mglobals.MINREADS
    min_npos_x_exp = mglobals.MINPOS
    in_data_filter = False
    cover = junc.coverage.toarray()[exp_idx]
    if junc.get_read_num(exp_idx) > min_read_x_exp and np.count_nonzero(cover) >= min_npos_x_exp:
        in_data_filter = True
    return in_data_filter


def lsv_detection(gene_list, chrom, logging=None):

    num_ss_var = [[0]*20, [0]*20, 0]

    const_set = [set() for xx in range(mglobals.num_experiments)]
    lsv_list = [[] for xx in range(mglobals.num_experiments)]
    jun = {}
    for xx in mglobals.tissue_repl.keys():
        jun[xx] = set()

    for strand, glist in gene_list.items():
        for gn in glist:

            gn.check_exons()
            count = gn.get_read_count().sum()
            if count == 0:
                continue
            mat, exon_list, tlb, var_ss = gn.get_rnaseq_mat(const_set, lsv=True)
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
            #utils.print_junc_matrices(mat, tlb=tlb, fp=True)
            SS, ST = lsv_matrix_detection(mat, tlb, (False, False, False), vip)
            dummy = {}
            for name, ind_list in mglobals.tissue_repl.items():
                dummy[name] = [[], []]

            for lsv_index, lsv_lst in enumerate((SS, ST)):
                lsv_type = (SSOURCE, STARGET)[lsv_index]
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
                        counter = 0
                        e_data = 0
                        for jj in jlist:
                            for exp_idx in ind_list:
                                if reliable_in_data(jj, exp_idx):
                                    counter += 1
                            if counter < 0.1*len(ind_list):
                                continue
                            e_data += 1
                            jun[name].add(jj)
                        if e_data == 0:
                            continue

                        dummy[name][lsv_index].append(lsv_in)
                        # for exp_idx in ind_list:
                        #     for lsvinlist in lsv_list[exp_idx]:
                        #         if lsv_in.is_equivalent(lsvinlist):
                        #             break
                        #     else:
                        #         if lsv_in.get_junctions_list().shape[0] >= 2:
                        #             lsv_list[exp_idx].append(lsv_in)

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

        # cand = range(ii+1, len(exon_to_ss))
        # if ii > 0:
        #     pre_lsv = exon_to_ss[ii-1]
        #     for ex_idx, ex in enumerate(cand):
        #         pt = exon_to_ss[ex]
        #         junc_cand = mat[lsv[1][0]:lsv[1][-1]+1, pt[0][0]:pt[0][-1]+1]
        #         if np.count_nonzero(junc_cand) < 1:
        #             continue
        #         to_trgt = mat[: pre_lsv[1][0]+1, pt[0][0]:pt[0][-1]+1]
        #         if np.count_nonzero(to_trgt) > 0 and not ex in vip_set:
        #             ss_valid = False
        #             break

        if ss_valid and np.count_nonzero(ss) > 1:
            lsv_list[0].append(ii)

    for ii in range(1, len(exon_to_ss)):
        lsv = exon_to_ss[ii]
        #Single Targe detection
        st = mat[:, lsv[0][0]:lsv[0][-1]+1]
        st_valid = True

        # cand = range(0, ii)
        # if ii+1 < len(exon_to_ss):
        #     post_lsv = exon_to_ss[ii+1]
        #     for ex_idx, ex in enumerate(cand):
        #         pt = exon_to_ss[ex]
        #         junc_cand = mat[pt[1][0]:pt[1][-1]+1, lsv[0][0]:lsv[0][-1]+1]
        #         if np.count_nonzero(junc_cand) < 1:
        #             continue
        #         from_src = mat[pt[1][0]:pt[1][-1]+1, post_lsv[0][0]:]
        #         if np.count_nonzero(from_src) > 0 and not ex in vip_set:
        #             st_valid = False
        #             break

        if st_valid and np.count_nonzero(st) > 1:
            lsv_list[1].append(ii)

    return lsv_list
