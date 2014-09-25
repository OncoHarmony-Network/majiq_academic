import numpy as np
from utils import utils
from lsv import SSOURCE, STARGET
import mglobals


def __junction_filter_check(junc):

    ''' Check total_read >= 10 '''
    junc_filter = False
#    if junc.readN[exp_idx] < 10 : return False 
    count = 0
    dummy = junc.coverage.toarray()
    for exp_idx in range(mglobals.num_experiments):
        cov_cnt = 0
        if junc.readN[exp_idx] >= 10:
            for val in dummy[exp_idx]:
                if val > 0:
                    cov_cnt += 1
            if cov_cnt < 3:
                continue
            count += 1
        else:
            continue
        if count > (0.1 * mglobals.num_experiments):
            junc_filter = True
            break
#    print "JUNC:", count
    return junc_filter


def __reliable_in_data(junc, exp_idx):
    
    min_read_x_exp = mglobals.MINREADS
    min_npos_x_exp = mglobals.MINPOS
    in_data_filter = False
    cover = junc.coverage.toarray()[exp_idx]
    if junc.readN[exp_idx] > min_read_x_exp and np.count_nonzero(cover) >= min_npos_x_exp:
        in_data_filter = True
    return in_data_filter


def __total_ss_minreads(junc_mat, minreads=5):
#    print " MIN READS on "
#    print junc_mat

    js = set()
    for jlst in junc_mat:
        for jj in jlst:
            if jj is None:
                continue
            if jj.readN.sum() > minreads:
                js.add(jj)

    return len(js)


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
                                if __reliable_in_data(jj, exp_idx):
                                    counter += 1
                            if counter < 0.1*len(ind_list):
                                continue
                            e_data += 1
                            jun[name].add(jj)
                        if e_data == 0:
                            continue
                        for exp_idx in ind_list:
                            for lsvinlist in lsv_list[exp_idx]:
                                if lsv_in.is_equivalent(lsvinlist):
                                    break
                            else:
                                if lsv_in.get_junctions_list().shape[0] >= 2:
                                    lsv_list[exp_idx].append(lsv_in)

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
    for ii in range(1, len(exon_to_ss) - 1):
        lsv = exon_to_ss[ii]
        pre_lsv = exon_to_ss[ii-1]
        post_lsv = exon_to_ss[ii+1]
        #Single Source detection
        ss = mat[lsv[1][0]:lsv[1][-1]+1, :]
        ss_valid = True
        cand = range(ii+1, len(exon_to_ss))
        for ex_idx, ex in enumerate(cand):
            pt = exon_to_ss[ex_idx]
            junc_cand = mat[lsv[1][0]:lsv[1][-1]+1, pt[0][0]:pt[0][-1]+1]
            if np.count_nonzero(junc_cand) < 1:
                continue
            to_trgt = mat[: pre_lsv[1][0]+1, pt[0][0]:pt[0][-1]+1]
            if np.count_nonzero(to_trgt) > 0 and not ex_idx in vip_set:
                ss_valid = False
                break

        if ss_valid and np.count_nonzero(ss) > 1:
            lsv_list[0].append(ii)

        #Single Targe detection
        st = mat[:, lsv[0][0]:lsv[0][-1]+1]
        st_valid = True
        cand = range(0, ii)
        for ex_idx, ex in enumerate(cand):
            pt = exon_to_ss[ex_idx]
            junc_cand = mat[pt[1][0]:pt[1][-1]+1, lsv[0][0]:lsv[0][-1]+1]
            if np.count_nonzero(junc_cand) < 1:
                continue
            from_src = mat[pt[1][0]:pt[1][-1]+1, post_lsv[0][0]:]
            if np.count_nonzero(from_src) > 0 and not ex_idx in vip_set:
                st_valid = False
                break

        if st_valid and np.count_nonzero(st) > 1:
            lsv_list[1].append(ii)

    return lsv_list