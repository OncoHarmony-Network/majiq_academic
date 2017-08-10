import pdb
from voila.tools import Tool
from voila.tools.find_binary_lsvs import find_binary_lsv_ids
from voila.tools.utils import io_caleb
from voila.tools.utils.most_common_elem import most_common
from voila.tools.utils.percent_through_list import percent_through_list
from voila.utils.voila_log import voila_log
import numpy as np
import copy
import pickle as pkl
import os

# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


LOG = voila_log()


class ThisisNonRedundantSets(Tool):
    help = 'Identify how many non-redundant AS events there are, whereby LSVs that are ' \
           'connected by \'utilized junctions\' are considered part of the same splicing event.' \
           '\'Utilized junctions\' means dPSI >= some threshold. (And/or Prob(dPSI) and PSI if you want).'

    def arguments(self):
        parser = self.get_parser()
        mutually_excl_grp = parser.add_mutually_exclusive_group(required=True)
        mutually_excl_grp.add_argument('--directory',
                                       '-d',
                                       type=str,
                                       help='Directory where voila texts are. (searches recursively)')
        mutually_excl_grp.add_argument('--file_list',
                                       '-fl',
                                       type=str,
                                       help='A file with a list of voila text file paths (one line each).')
        help_mes = "Output to save file to (include full path)"
        parser.add_argument('--outfile',
                            type=str,
                            help=help_mes)
        help_mes = "E(dPSI) threshold by which to call junctions as utilized"
        parser.add_argument('--dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.1)
        # TODO if requested...
        # help_mes = "Prob(dPSI) threshold by which to call junctions as utilized"
        # parser.add_argument('--prob_dpsi_thresh',
        #                     type=float,
        #                     help=help_mes,
        #                     default=0.0)
        help_mes = "E(PSI) threshold by which to call junctions as utilized. 1.0 essentially means" \
                   "don't use PSI..."
        parser.add_argument('--psi_thresh',
                            type=float,
                            help=help_mes,
                            default=1.0)
        help_mes = "Flag: don't consider intron retention LSVs"
        parser.add_argument('--no_ir',
                            action='store_true',
                            help=help_mes,
                            default=False)
        help_mes = "Flag: save as pickle file instead of list"
        parser.add_argument('--ferment',
                            action='store_true',
                            help=help_mes,
                            default=False)
        help_mes = 'Optional pattern matching to identify the voila text files'
        parser.add_argument('-p',
                            '--pattern',
                            default="*tsv",
                            type=str,
                            help=help_mes)
        return parser

    def run(self, args):
        # this is for code readability, not efficiency
        consider_ir = True
        if args.no_ir:
            consider_ir = False
        imported = io_caleb.quick_import(input=args.directory if args.directory else args.file_list,
                                         cutoff_d_psi=0,
                                         cutoff_prob=0,
                                         pattern=args.pattern,
                                         keep_ir=consider_ir)
        non_red_res = non_redundant_set(imported,
                                        cutoff_dpsi=args.dpsi_thresh,
                                        cutoff_psi=args.psi_thresh)
        if args.outfile:
            out_path = os.path.abspath(args.outfile)
            if args.ferment:
                LOG.info("Writing non-redundant set data as pickle file to %s" % out_path)
                pkl.dump(non_red_res, open(out_path, "wb"))
            else:
                LOG.info("Writing non-redundant sets line-by-line to %s" % out_path)
                non_red_list = non_red_sets_to_list(non_red_res)
                with open(out_path, "w") as handle:
                    for line in non_red_list:
                        line.sort()
                        print("|".join(line), file=handle)
        else:
            non_red_list = non_red_sets_to_list(non_red_res)
            for line in non_red_list:
                line.sort()
                print("|".join(line))


def non_redundant_set(data,
                      cutoff_dpsi=0.1,
                      cutoff_psi=1.0,
                      save_blanked_structure=False,
                      return_numdpsis_dat=False):
    """
    Identify how many non-redundant AS events there are, whereby LSVs
        that are connected by 'utilized junctions' are considered part
        of the same splicing event.

            "Utilized junctions": utilized meaning, dPSI/PSI cutoff
                met in at least one comparison/group in the data.

        Arguments:
            data : quick import...
            cutoff_dpsi: how much must the junction dPSI be (in at least
                one comparison) in order to be considered utilized?
            cutoff_psi: If PSI of junciton is greater than this cutoff in
                at least one sample, it is considered utilized?
                Default is 1 : no junction makes the PSI cutoff >1.0
            save_blanked_structure: if True, don't revert the Data.
            However, need to also return the blanked_lsvs_dict so you
            know which LSVs were blanked.
            return_numdpsis_dat: boolean

    """
    # First need to fill in LSV Gaps. deltapsi comparisons don't overlap 100%
    # with respect to LSVs. This will add 'blank' (dPSI, PSI, Prob =0) LSVs
    # to the Data. I'm doing it InPlace, so it alters the inut object.
    # No worries, though, because I'll revert the Data to its original state.
    blanked_lsvs_dict = io_caleb.impute_missing_lsvs(data=data, in_place=True, warnings=False)
    LOG.info("Finished filling in the gaps, running non-redundant algorithm...")
    if return_numdpsis_dat:
        nr_connected_lsvs, nr_numdpsis = get_connected_lsvs_by_junc(data=data,
                                                                    Cutoff_dPSI=cutoff_dpsi,
                                                                    Cutoff_PSI=cutoff_psi,
                                                                    ret_numpdsis_data=return_numdpsis_dat)
    else:
        nr_connected_lsvs = get_connected_lsvs_by_junc(data=data,
                                                       Cutoff_dPSI=cutoff_dpsi,
                                                       Cutoff_PSI=cutoff_psi,
                                                       ret_numpdsis_data=return_numdpsis_dat)
    LOG.info("Accounting for LSVs from different genes that overlap according to genomic coordinates ...")
    by_type = most_lsvs_same_gene(connected_lsv_s=nr_connected_lsvs)
    if not save_blanked_structure:
        LOG.info("Reverting Data to original state....")
        for comparison in data.keys():
            for lsv_id in blanked_lsvs_dict[comparison]:
                data[comparison].pop(lsv_id)
    LOG.info("Done!!!\n\n\n")
    s = len(by_type["singles"])
    ns = len(by_type["singles_all_lsvs"])
    d = len(by_type["twos"])
    nd = len(by_type["twos_lsvs"])
    tp = len(by_type["three_plus"])
    ntp = len(by_type["three_plus_lsvs"])
    total = s + d + tp
    total_lsvs = ns + nd + ntp
    summary = "%s total LSVs met junction utilization cutoffs...\n" % total_lsvs
    summary += "%s total non-redundant splicing sets\n" % total
    summary += "%s singles (%s LSVs), %s doubles (%s LSVs), %s three_pluses (%s LSVs) ...\n" % (s,
                                                                                            ns,
                                                                                            d,
                                                                                            nd,
                                                                                            tp,
                                                                                            ntp)
    LOG.info(summary)
    by_type["summary"] = summary
    if return_numdpsis_dat:
        if save_blanked_structure:
            return by_type, blanked_lsvs_dict, nr_numdpsis
        return by_type, nr_numdpsis
    else:
        if save_blanked_structure:
            return by_type, blanked_lsvs_dict
        return by_type


def non_red_sets_to_list(non_red_data):
    """ Convert dictionary results from non_red_data to
        a list of non-red LSV sets.

    :param non_red_data:
    :return: [[non red set 1], [non red set 2], [etc] ]
    """
    singles = list(non_red_data["singles"].values())
    twos = list(non_red_data["twos"].values())
    three_plus = list(non_red_data["three_plus"].values())
    all_non_red_sets = list()
    all_non_red_sets.extend(singles)
    all_non_red_sets.extend(twos)
    all_non_red_sets.extend(three_plus)
    return all_non_red_sets


def most_lsvs_same_gene(connected_lsv_s,
                        return_shared_junc_genes=False):
    """
    Some lsvs that are connected by a junction turn out to
        to be due to overlapping genes. Thus, this function teases
        out how which connected LSVs correspond to connections in
        the same gene. Removes LSVs that are redundant.

        Return_shared_junc_genes: if True, return list of junctions that
            multiple LSVs (in more than 1 gene) share in the Connected
            LSV network.
    """
    these_comps = list(connected_lsv_s.keys())
    lengths = [len(connected_lsv_s[x]) for x in these_comps]
    lengths = np.array(lengths)
    singles = np.where(lengths == 1)[0]
    singles = singles.tolist()
    singles_lsvs = list()
    for single in singles:
        singles_lsvs.extend(connected_lsv_s[these_comps[single]])
    results = dict()
    twos = list()
    two_lsvs = list()
    three_plus = list()
    three_lsvs = list()
    two_plus_connected = np.where(lengths > 1)[0]
    two_plus_connected = two_plus_connected.tolist()
    genes_share_junc = list()
    for two_plus in two_plus_connected:
        # pdb.set_trace()
        connected = connected_lsv_s[these_comps[two_plus]]
        gene_ids = [x.split(":")[0] for x in connected]
        if len(set(gene_ids)) > 1:
            genes_share_junc.append(these_comps[two_plus])
        # gene_with_mult_lsvs = [x for x in gene_ids if gene_ids.count(x) > 1]
        gene_with_max_lsvs = most_common(gene_ids)
        n_max_lsvs_one_gene = gene_ids.count(gene_with_max_lsvs)
        # if count_of_most_common>1:
        #     gene_with_max_lsvs = max([x for x in gene_ids if gene_ids.count(x) > 1])
        if n_max_lsvs_one_gene == 1:
            singles.append(two_plus)  # there is only 1 lsv from each gene in this connection
            singles_lsvs.extend(connected)
            continue
        elif n_max_lsvs_one_gene == 0:
            raise RuntimeError("No idea what happened...")
        n_max_lsvs_one_gene = gene_ids.count(gene_with_max_lsvs)
        if n_max_lsvs_one_gene > 2:
            # print "==="
            # print connected
            # print two_plus
            three_plus.append(two_plus)
            three_lsvs.extend(connected)
            # print "==="
        if n_max_lsvs_one_gene == 2:
            twos.append(two_plus)
            two_lsvs.extend(connected)
    results["singles"] = {these_comps[k]: connected_lsv_s[these_comps[k]] for k in singles}
    results["singles_all_lsvs"] = singles_lsvs
    results["twos"] = {these_comps[k]: connected_lsv_s[these_comps[k]] for k in twos}
    results["twos_lsvs"] = two_lsvs
    results["three_plus"] = {these_comps[k]: connected_lsv_s[these_comps[k]] for k in three_plus}
    results["three_plus_lsvs"] = three_lsvs
    if return_shared_junc_genes:
        dict_res = dict()
        dict_res["types"] = results
        dict_res["genes_share_junc"] = genes_share_junc
        return dict_res
    return results



def get_connected_lsvs_by_junc(data,
                               Cutoff_dPSI=0.1,
                               Cutoff_PSI=1.0,
                               ret_numpdsis_data=False):
    """
        Given LSV quick import, identify which LSVs are connected
            using utilized junctions as edges and LSVs as nodes.

            Utilized: may use dPSI OR PSI to define whether a junction
            is utilized.

            Cutoff_dPSI: how much must the junction dPSI be (in at least
                one comparison) in order to be considered utilized?
            Cutoff_PSI: What PSI must the junciton be greater than in order
                to be considered utilized? If set to 1, then no junction
                makes the PSI cutoff.
            ret_numpdsis_data: boolean

    """
    # sig_juncs = get_sorted_juncs(Data)
    master_junc_dict = dict()
    all_lsvs = io_caleb.get_shared_lsv_ids(data)
    junc_lsv_dicts = get_junc_lsv_dicts(data,
                                        cutoff_dpsi=Cutoff_dPSI,
                                        cutoff_psi=Cutoff_PSI,
                                        return_numdpsis=ret_numpdsis_data)
    junc_to_lsv = junc_lsv_dicts["junc_lsv_dict"]
    lsv_to_junc = junc_lsv_dicts["lsv_junc_dict"]
    sig_juncs = list(junc_to_lsv.keys())
    LOG.info("Recursively processing %s sig junctions to build connectivity graph ..." % (len(sig_juncs)))
    indeces_at_x_percent = percent_through_list(len(sig_juncs), 0.01)
    i = 0.0
    n_to_start = len(sig_juncs)
    while len(sig_juncs) > 0:
        if i > 0.0 and i in indeces_at_x_percent:
            LOG.info(str(indeces_at_x_percent[i]) + "% of sig juncs processed...")
        next_junc = sig_juncs.pop()
        conn_lsvs = connected_lsvs(all_lsvs, sig_juncs, junc_to_lsv, lsv_to_junc, next_junc)
        master_junc_dict[next_junc] = conn_lsvs
        i = n_to_start - len(sig_juncs)
    if ret_numpdsis_data:
        return master_junc_dict, junc_lsv_dicts["numdPSIs_data"]
    return master_junc_dict


def connected_lsvs(all_lsvs, sig_juncs, junc_to_lsv, lsv_to_Junc, start_junc):
    """
    Return all the LSVs connected to the start junc, recursively.
    """
    conn_lsvs = junc_to_lsv[start_junc]
    associated_juncs = list()
    associated_lsvs = list()
    # if Start_junc == 'chr11_+_8705628-8706265':
    #     pdb.set_trace()
    # if Start_junc == 'chr11_+_8705628-8706370':
    #     pdb.set_trace()
    for found_lsv in conn_lsvs:
        if found_lsv not in all_lsvs:
            continue
        else:
            associated_lsvs.append(found_lsv)
            all_lsvs.remove(found_lsv)
            found_lsv_juncs = lsv_to_Junc[found_lsv]
            found_lsv_juncs = found_lsv_juncs.tolist()
            associated_juncs.extend(found_lsv_juncs)
    associated_juncs = [x for x in associated_juncs if x in sig_juncs]  # only keep unseen juncs
    for found_junc in associated_juncs:
        if found_junc not in sig_juncs:
            continue
        else:
            sig_juncs.remove(found_junc)
            recursive_call = connected_lsvs(all_lsvs,
                                            sig_juncs,
                                            junc_to_lsv,
                                            lsv_to_Junc,
                                            found_junc)
            associated_lsvs.extend(recursive_call)
    return associated_lsvs


def get_sorted_juncs(Data):
    """
    Return list of junctions sorted by the maximum
        dPSI they have across all comparisons in
        the Data.
    """
    master_junc_weight_dict = dict()
    if io_caleb.check_is_lsv_dict(Data, da_bool=True):
        master_junc_weight_dict = get_junc_weights(Data)
    elif io_caleb.check_is_quick_import(Data):
        shared_ids = io_caleb.get_shared_lsv_ids(Data)
        comp_to_junc_weights = dict()
        all_comparisons = Data.keys()
        for comparison in all_comparisons:
            share_dict = io_caleb.lsv_dict_subset(Data[comparison], shared_ids)
            junc_to_weights = get_junc_weights(share_dict, Weights="dPSI")
            comp_to_junc_weights[comparison] = junc_to_weights
        all_juncs = comp_to_junc_weights[all_comparisons[0]].keys()

        for junc_coord in all_juncs:
            max_weight = max([comp_to_junc_weights[comp][junc_coord] for comp in all_comparisons])
            master_junc_weight_dict[junc_coord] = max_weight
    else:
        raise RuntimeError("Data needs to be LSV dict or quick import...")
    junc_coords = np.array(master_junc_weight_dict.keys())
    weights = np.array(master_junc_weight_dict.values())
    sort_inds_by_weight = np.argsort(weights)
    sort_inds_by_weight = sort_inds_by_weight[::-1]  # reverse order
    sorted_coords = junc_coords[sort_inds_by_weight]
    sorted_coords = sorted_coords.tolist()
    return sorted_coords


def get_master_junc_weights(Data, Weights="dPSI"):
    """
    Given quick import, look through all comparisons
        for the max Weight associated with a junc and
        return dictionary of junctions
        pointing at their maximum abs(dPSI) or PSI value.

        chr_strand_coords : Weight

        Arguments:
        Weight = "dPSI" or "PSI" <- which value to use?
    """
    io_caleb.check_is_quick_import(Data)
    comparisons = Data.keys()
    shared_ids = io_caleb.get_shared_lsv_ids(Data)
    shared_dicts = dict()
    for comparison in comparisons:
        share_dict = io_caleb.lsv_dict_subset(Data[comparison], shared_ids)
        junc_to_weights = get_junc_weights(share_dict, Weights=Weights)
        shared_dicts[comparison] = junc_to_weights
    all_juncs = shared_dicts[comparisons[0]].keys()
    max_dict = dict()
    for junc_key in all_juncs:
        max_weight = max([shared_dicts[comp][junc_key] for comp in comparisons])
        max_dict[junc_key] = max_weight
    return max_dict


def get_junc_weights(Data, Weights="dPSI"):
    """
    Given LSV dict, return dictionary of junctions
        pointing at their maximum abs(dPSI) or PSI value.

        Arguments:
        Weight = "dPSI" or "PSI" <- which value to use?

        PS: if dPSI, uses abs(dPSI)
    """
    if Weights != "dPSI" and Weights != "PSI":
        raise RuntimeError("Weights mst be 'dPSI' or 'PSI'")
    io_caleb.check_is_lsv_dict(Data)
    lsv_ids = io_caleb.get_lsv_ids(Data)
    junc_weight_dict = dict()
    all_juncs = dict()
    for lsv_id in lsv_ids:
        lsv = Data[lsv_id]
        chrm = lsv["chr"]
        strand = lsv["strand"]
        junc_coords = copy.copy(lsv["Junctions coords"])
        if Weights == "dPSI":
            weights = copy.copy(lsv["E(dPSI) per LSV junction"])
        else:
            weights_1 = copy.copy(lsv["E(PSI)1"])
            weights_2 = copy.copy(lsv["E(PSI)2"])
            weights = list()
            for weight_1, weight_2 in zip(weights_1, weights_2):
                weights.append(max([weight_1, weight_2]))
        new_coords = list()
        for coord in junc_coords:
            new_coords.append(chrm + "_" + strand + "_" + coord)
        for weight, coord in zip(weights, new_coords):
            if coord in all_juncs:
                all_juncs[coord] = max(all_juncs[coord], abs(weight))
            else:
                all_juncs[coord] = abs(weight)
    return all_juncs


def get_junc_lsv_dicts(data,
                       cutoff_dpsi=0.1,
                       cutoff_psi=0.05,
                       return_numdpsis=False):
    """
        Return a dictionary of junctions pointing at lists of LSV IDs
            that use them significantly. Also return dict of LSV IDs
            pointing at all their junctions (utilized only).

    """
    io_caleb.check_is_quick_import(data)
    LOG.info("Getting all dPSI data ...")
    this_num_dpsi = io_caleb.get_num_d_psi(data)
    LOG.info("Getting all PSI data ...")
    this_nump_psi = io_caleb.get_num_psi(data)
    LOG.info("Getting all junction coordinates ...")
    lsv_junc_dict = get_junc_coords(data)
    LOG.info("Running find_binary_LSV_IDs() ...")
    sig_junc_dict = find_binary_lsv_ids(num_d_psi=this_num_dpsi,
                                        num_psi=this_nump_psi,
                                        cutoff_d_psi=cutoff_dpsi,
                                        cutoff_psi=cutoff_psi,
                                        return_sig_juncs=True)
    junc_lsv_dict = dict()
    for lsv_id in sig_junc_dict.keys():
        junc_loci = lsv_junc_dict[lsv_id]
        try:
            sig_junc_bools = sig_junc_dict[lsv_id] > 0
            sig_coords = junc_loci[sig_junc_bools]
            lsv_junc_dict[lsv_id] = sig_coords
            if len(lsv_junc_dict[lsv_id]) == 0:
                del lsv_junc_dict[lsv_id]
        except:
            pdb.set_trace()
        for sig_coord in sig_coords:
            if sig_coord in junc_lsv_dict:
                junc_lsv_dict[sig_coord].append(lsv_id)
            else:
                junc_lsv_dict[sig_coord] = list()
                junc_lsv_dict[sig_coord].append(lsv_id)
    res = dict()
    res["junc_lsv_dict"] = junc_lsv_dict
    res["lsv_junc_dict"] = lsv_junc_dict
    if return_numdpsis:
        res["numdPSIs_data"] = this_num_dpsi
    return res


def get_junc_coords(Data):
    """
    Given LSV dict info (quick import, or LSV dict),
     return dict of lsv_ids pointing at array of junction loci:
        'chr_strand_start_end'
    """
    if io_caleb.check_is_quick_import(Data, the_bool=True):
        junc_dict = dict()
        for comp in Data:
            junc_dict.update(get_junc_coords(Data[comp]))
        return junc_dict
    io_caleb.check_is_lsv_dict(Data)
    junc_dict = dict()
    lsv_ids = io_caleb.get_lsv_ids(Data)
    for lsv_id in lsv_ids:
        lsv = Data[lsv_id]
        chrm = lsv["chr"]
        strand = lsv["strand"]
        junc_coords = copy.copy(lsv["Junctions coords"])
        new_coords = list()
        for coord in junc_coords:
            new_coords.append(chrm + "_" + strand + "_" + coord)
        new_coords = np.array(new_coords)
        junc_dict[lsv_id] = new_coords
    return junc_dict
