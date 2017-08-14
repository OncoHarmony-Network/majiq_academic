from voila.tools import Tool
from voila.tools import io_voila_caleb
import os
import itertools
import math
import numpy as np
import pandas as pa
import copy
import operator
import pdb


class ThisisFindBinaryLSVs(Tool):
    help = 'Given a list of voila dPSI txt files, return LSVs that are binary-like'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory where voila texts are.')
        help_mes = "dPSI threshold by which to call junctions as changing"
        parser.add_argument('--dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.1)
        help_mes = "Prob(dPSI) threshold by which to call junctions as changing"
        parser.add_argument('--prob_dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.0)
        help_mes = "Flag: don't consider IR LSVs"
        parser.add_argument('--no_ir',
                            action='store_true',
                            help=help_mes,
                            default=False)
        help_mes = 'Optional pattern matching to identify the voila text files'
        parser.add_argument('-p',
                            '--pattern',
                            default="*tsv",
                            type=str,
                            help=help_mes)
        help_mes = 'Method: \'count\' or \'sum_to_95\''
        parser.add_argument('-m',
                            '--method',
                            choices={"sum_to_95", "count"},
                            default="sum_to_95",
                            type=str,
                            help=help_mes)
        help_mes = "Flag: just return the LSV IDs?"
        parser.add_argument('--just_ids',
                            action='store_true',
                            help=help_mes,
                            default=False)
        help_mes = "Flag: If you don't care whether the binary-like LSVs are reciprocating"
        parser.add_argument('--non_reciprocating',
                            action='store_false',
                            help=help_mes,
                            default=True)
        help_mes = "Flag: If you don't want to impute missing values with 0"
        parser.add_argument('--dont_impute',
                            action='store_true',
                            help=help_mes,
                            default=False)
        return parser

    def run(self, args):
        # this is for code readability, not efficiency
        consider_ir = True
        if args.no_ir:
            consider_ir = False
        impute = True
        if args.dont_impute:
            impute = False
        imported = io_voila_caleb.quick_import(dir=args.directory,
                                               cutoff_d_psi=0,
                                               cutoff_prob=0,
                                               pattern=args.pattern,
                                               keep_ir=consider_ir)
        if impute:
            imported = io_voila_caleb.impute_missing_lsvs(data=imported,
                                                          impute_with=0)
        results = get_binary_lsvs(data=imported,
                                  method=args.method,
                                  cutoff_d_psi=args.dpsi_thresh,
                                  just_lsv_ids=args.just_ids,
                                  must_reciprocate=args.non_reciprocating)
        if args.prob_dpsi_thresh:
            results = io_voila_caleb.subset_significant(results,
                                                        cutoff_dpsi=0,
                                                        cutoff_prob=args.prob_dpsi_thresh,
                                                        keep_introns=consider_ir)
        if args.just_ids:
            print(results)


def get_binary_lsvs(data,
                    cutoff_d_psi=0.1,
                    cutoff_psi=1,
                    make_copy=False,
                    just_lsv_ids=False,
                    method="sum_to_95",
                    must_reciprocate=True):
    """
    Given a dictionary of LSV dictionaries, return a subset of the
        dictionary whereby the LSVs all look functionally binary.

    See find_binary_LSV_IDs for details.

    Arguments:
        data: quick import
        cutoff_d_psi: dpsi
        cutoff_psi: psi
        make_copy: if True, don't modify the original Data, make a copy of
            the subset. If False, the original Data will be overwritten.
        just_lsv_ids: if False, return LSV dicts
        must_reciprocate: Trye/False, binary-like LSVs must have a junc go up AND down
        method: 'sum_to_95' or 'count' ... see find_binary_lsv_ids and find_binary_lsvs_95
    """
    io_voila_caleb.check_is_quick_import(data)

    print("Getting num_d_psi data ...")

    num_d_psi_data = get_num_d_psi(data)
    if cutoff_psi < 1:
        print("Getting num_psi data ...")
        numPSIs_data = get_num_psi(data)
    else:
        numPSIs_data = "doesnt matter"
    print("Categorizing LSVs from the following comparisons ...")
    print(list(data.keys()))
    if method == "count":
        results = find_binary_lsv_ids(num_d_psi_data,
                                      numPSIs_data,
                                      cutoff_d_psi=cutoff_d_psi,
                                      cutoff_psi=cutoff_psi,
                                      return_bi_iis=True,
                                      must_reciprocate=must_reciprocate)
    elif method == "sum_to_95":
        results = find_binary_lsvs_95(num_d_psi_data,
                                      numPSIs_data,
                                      threshold=cutoff_d_psi,
                                      by_d_psi=True,
                                      by_psi=False,
                                      return_bi_iis=True,
                                      must_reciprocate=must_reciprocate)
    else:
        raise ValueError("Expected 'count' or 'sum_to_95', not %s" % method)
    binary_ids = results["binary"]

    if just_lsv_ids:
        return binary_ids
    binary_indices = results["binary_indices"]

    if make_copy:
        return_dict = dict()
    for comparison in list(data.keys()):
        if make_copy:
            return_dict[comparison] = io_voila_caleb.lsv_dict_subset(dictionary=data[comparison],
                                                                     keys=binary_ids,
                                                                     save_LSV_data=True,
                                                                     new_sub_key="binary_indices",
                                                                     new_values=binary_indices)
        else:
            data[comparison] = io_voila_caleb.lsv_dict_subset(dictionary=data[comparison],
                                                              keys=binary_ids,
                                                              save_LSV_data=True,
                                                              new_sub_key="binary_indices",
                                                              new_values=binary_indices)
    if make_copy:
        return return_dict
    else:
        return data


def find_binary_lsv_ids(num_d_psi,
                        num_psi,
                        cutoff_d_psi=0.1,
                        cutoff_psi=1,
                        return_other_ids=False,
                        return_bi_iis=False,
                        return_all_iis=False,
                        return_sig_juncs=False,
                        debug=False,
                        must_reciprocate=True):
    """
    Looking only at junctions that meet Cutoff_PSI and Cutoff_dPSI,
        identify which LSVs have 0, 1, 2, or 3+ junctions changing
        across all conditions.

    Cutoff_PSI is used to make sure that LSVs that have non-changing
        junctions (low dPSI), but are still included at Cutoff_PSI amount
        are considered "used" when determining binary-like behavior of the LSV.
        Set this to 1 if you only want to categorize LSVs by dPSI.

    Arguments
        Cutoff_dPSI: default dPSI cutoff is 0.1.
        Cutoff_PSI: defualt PSI cutoff is 1
        Return_other_ids: if True, return LSV IDs that have:
                     ==2, >2 , ==1, ==0
            junctions changing accross all comparisons.
        Return_binary_indices: if True, also return a dictionary
            of binary lsv_ids pointing at the indices for each
            junction that is functionally binary
        Return_all_indeces: if True, return dictionary of lsv_ids
            pointing at the indeces showing Cutoff_dPSI in any comparison
        Return_sig_juncs: if True, return dict[lsv_ids] pointing at
            arrays indicating how many junctions found to be 'significant'
            Note: for now, this overrides all other things returned ... TODO: fixthatshit
        debug: boolean, use set trace?
        must_reciprocate: Bool: do binary-like LSVs need a junc going up and the other going down?

        If the #True column has exactly two junctions > 0,
        then the LSV is considered functionally binary.

                |A vs B|A vs C  | ... |
    Junction:   | dPSI | dPSI   | ... |
            0   |>0.1? |>0.1?   | ... | #True
            1   |>0.1? |>0.1?   | ... | #True
            2   |>0.1? |>0.1?   | ... | #True
            ... |   ...|  ...   | ... | ...
    """
    lsv_ids = list(num_d_psi.keys())
    binary_ids = list()  # exacly 2 junctions over cutoff
    binary_indices = list()
    complex_over_ids = list()  # >2 junctions
    complex_single_ids = list()  # single junction
    zero_over_ids = list()  # zero junctions
    all_indices = dict()
    sig_juncs_dict = dict()
    print("Counting how many juncs utilized in %s LSVs ..." % len(lsv_ids))
    i = 1.0
    indeces_at_10_percent = percent_through_list(lsv_ids, 0.1)
    for lsv_id in lsv_ids:
        if i > 0.0 and i in indeces_at_10_percent:
            print(str(indeces_at_10_percent[i]) + "% of juncs looked at...")
        i += 1.0
        this_num_d_psi = num_d_psi[lsv_id]

        # identify which junctions have PSI > Cutoff_dPSI
        dpsi_over_cutoff = abs(this_num_d_psi) > cutoff_d_psi
        if not isinstance(num_psi, str):
            this_num_psi = num_psi[lsv_id]

            # identify which junctions have PSI > Cutoff_PSI
            psi_over_cutoff = np.sum(this_num_psi > cutoff_psi, axis=1)
            psi_over_cutoff = np.zeros(this_num_d_psi.shape).T + psi_over_cutoff
            psi_over_cutoff = psi_over_cutoff.T > 0
            # make sure all junctions with dPSI and PSI over the cutoff
            # are 'True'
            dpsi_over_cutoff = psi_over_cutoff + dpsi_over_cutoff
        # axis=1, meaning sum # of Trues in each row (each junction)
        sum_truths = np.sum(dpsi_over_cutoff, axis=1)
        sig_juncs_dict[lsv_id] = sum_truths
        number_junctions_over_cutoff = sum(sum_truths > 0)
        if number_junctions_over_cutoff == 2:
            if debug:
                pdb.set_trace()
            doesnt_reciprocate = False
            if must_reciprocate:
                abs_junc_maxes = np.max(abs(this_num_d_psi), axis=1)
                # reverse sorted indices
                top_indices = np.argpartition(abs_junc_maxes, -2)[-2:]
                if not does_reciprocate(num_d_psi=this_num_d_psi, top_indices=top_indices):
                    doesnt_reciprocate = True
            if doesnt_reciprocate:
                complex_over_ids.append(lsv_id)
            else:
                binary_ids.append(lsv_id)
                if return_bi_iis:
                    # which rows have a total of 0 True?
                    indices = index_all(sum_truths.tolist(), 0, nottt=True)
                    binary_indices.append(indices)
        elif number_junctions_over_cutoff > 2:
            complex_over_ids.append(lsv_id)
        elif number_junctions_over_cutoff == 1:
            complex_single_ids.append(lsv_id)
        else:  # == 0
            zero_over_ids.append(lsv_id)
        all_indices[lsv_id] = index_all(sum_truths.tolist(), 0, nottt=True)
    if return_sig_juncs:
        return sig_juncs_dict
    results = dict()
    results["binary"] = binary_ids
    if return_other_ids:
        results["complex"] = complex_over_ids
        results["single"] = complex_single_ids
        results["non_sig"] = zero_over_ids
    if return_bi_iis:
        results["binary_indices"] = binary_indices
    if return_all_iis:
        results["all_indices"] = all_indices
    n_0 = len(zero_over_ids)
    n_1 = len(complex_single_ids)
    n_2 = len(binary_ids)
    n_comp = len(complex_over_ids)
    print("%s non-sig, %s single junc, %s binary-like, %s complex LSVs categorized." % (n_0, n_1, n_2, n_comp))
    return results


def find_binary_lsvs_95(num_d_psi,
                        num_psi,
                        threshold=0.1,
                        by_d_psi=True,
                        by_psi=False,
                        return_other_ids=False,
                        return_bi_iis=False,
                        return_all_iis=False,
                        # return_sig_juncs=False,
                        must_reciprocate=True):  # meaning, top 2  juncs are positive AND negative
    """
    if Sum(abs(top two junctions dPSIs)) / Sum(abs(dPSI)) >= 0.95, its binary.
        Use the maximum dPSI seen across all comparisons to compute above thing.TODO: ask Yoseph if that's cool.
    """
    lsv_ids = list(num_d_psi.keys())
    binary_ids = list()  # exacly 2 junctions over cutoff
    binary_indices = list()
    complex_over_ids = list()  # >2 junctions
    complex_single_ids = list()  # single junction
    zero_over_ids = list()  # zero junctions
    all_indices = dict()
    sig_juncs_dict = dict()
    print("%s shared LSVs being categorized..." % len(lsv_ids))
    for lsv_id in lsv_ids:
        # if lsv_id == 'ENSG00000003756:50131151-50131763:target':
        #     #pass
        #     pdb.set_trace()
        this_num_d_psi = num_d_psi[lsv_id]
        if not isinstance(num_psi, str):
            this_num_psi = num_psi[lsv_id]
        else:
            this_num_psi = "not using psi..."
        which_junctions_utilized = np.zeros(this_num_d_psi.shape[0], dtype=bool)
        top = 0  # initialize
        # Use PSI to determing binary-ness:
        if by_psi:
            print("Gotta be straight with you, Caleb was lazy and didn't implement "
                  "the thresh option for binary PSI... fix this")
            abs_junc_maxes = np.max(abs(this_num_psi), axis=1)  # max per row (junction)
            i = -1
            perc_of_total = 0
            while perc_of_total < 0.95 and abs(i) <= len(abs_junc_maxes):
                # print i
                # indices in reverse order of the max values
                top_indices = np.argpartition(abs_junc_maxes, i)[i:]
                top = abs_junc_maxes[top_indices]  # top values, rev order
                perc_of_total = np.sum(top) / np.sum(abs_junc_maxes)
                # print perc_of_total
                i -= 1
            if i == len(abs_junc_maxes) and perc_of_total < 0.95:
                pdb.set_trace()
            which_junctions_utilized[top_indices] = True
        if by_d_psi:
            # Use dPSI to determing binary-ness:
            abs_junc_maxes = np.max(abs(this_num_d_psi), axis=1)  # max per row (junction)
            i = -1
            perc_of_total = 0
            while perc_of_total < 0.95 and abs(i) <= len(abs_junc_maxes):
                # reverse sorted indices
                top_indices = np.argpartition(abs_junc_maxes, i)[i:]
                top = abs_junc_maxes[top_indices]  # top values
                perc_of_total = np.sum(top) / np.sum(abs_junc_maxes)
                # print perc_of_total
                i -= 1
            if i == len(abs_junc_maxes) and perc_of_total < 0.95:
                pdb.set_trace()
            which_junctions_utilized[top_indices] = True

        all_indices[lsv_id] = index_all(which_junctions_utilized.tolist(), 0, nottt=True)
        if not by_psi and not by_d_psi:
            raise RuntimeError("Hey, you need to pick PSI or dPSI...")

        # Check if any of the dPSI or PSI are meating your thresh
        over_thresh = True if max(top) >= threshold else False

        # axis=1, meaning sum # of Trues in each row (each junction)
        sum_truths = np.sum(which_junctions_utilized)
        if sum_truths == 2 and over_thresh:
            doesnt_reciprocate = False
            if must_reciprocate:
                if not does_reciprocate(num_d_psi=this_num_d_psi, top_indices=top_indices):
                    doesnt_reciprocate = True
            if doesnt_reciprocate:
                complex_over_ids.append(lsv_id)
            else:
                binary_ids.append(lsv_id)
                if return_bi_iis:
                    # which rows have a total of 0 True?
                    indices = index_all(which_junctions_utilized.tolist(), 0, nottt=True)
                    binary_indices.append(indices)
        elif sum_truths > 2 and over_thresh:
            complex_over_ids.append(lsv_id)
        elif sum_truths == 1 and over_thresh:
            complex_single_ids.append(lsv_id)
        elif over_thresh:  # ????
            pdb.set_trace()
        else:  # == 0
            zero_over_ids.append(lsv_id)

    results = dict()
    results["binary"] = binary_ids
    if return_other_ids:
        results["complex"] = complex_over_ids
        results["single"] = complex_single_ids
        results["non_sig"] = zero_over_ids
    if return_bi_iis:
        results["binary_indices"] = binary_indices
    if return_all_iis:
        results["all_indices"] = all_indices
    n_0 = len(zero_over_ids)
    n_1 = len(complex_single_ids)
    n_2 = len(binary_ids)
    n_comp = len(complex_over_ids)
    print("%s non-sig, %s single junc, %s binary-like, %s complex LSVs categorized." % (n_0,
                                                                                        n_1,
                                                                                        n_2,
                                                                                        n_comp))
    return results


def does_reciprocate(num_d_psi, top_indices):
    """

    :param num_d_psi: numpy array of dPSI (rows=junctions, cols=comparisons)
    :param top_indices: inidices of the biggest two abs(dPSI)
    :return: True or False
    """
    if min(np.max(num_d_psi, axis=1)[top_indices]) >= 0 or \
                    max(np.max(num_d_psi, axis=1)[top_indices]) <= 0:
        return False
    return True


def get_binary_type(lsv_a, lsv_b):
    """
    Returns "ME", "TC", or "Neither" for:
        Mutually exclusive,
        Tandem Cassette
        or
        Neither
    """
    io_voila_caleb.check_is_lsv(lsv_a)
    io_voila_caleb.check_is_lsv(lsv_b)
    if not lsv_a.has_key("binary_indices"):
        raise RuntimeError("Expected binary indices for A")
    if not lsv_b.has_key("binary_indices"):
        raise RuntimeError("Expected binary indices for B")
    # impossible to be mutually exclusive..
    if lsv_a["Reference_Type"] == lsv_b["Reference_Type"]:
        return "Neither"
    lsv_a_juncs = copy.copy(lsv_a['Junctions coords'])
    lsv_a_juncs = [lsv_a_juncs[x].split("-") for x in lsv_a["binary_indices"]]
    lsv_a_juncs = [map(int, x) for x in lsv_a_juncs]
    flattened_lsv_a_juncs = flatten_list(lsv_a_juncs)
    if len(set(flattened_lsv_a_juncs)) != 3:
        return "Neither"
    shared_junc_coord_lsv_a = list_duplicates(flattened_lsv_a_juncs)[0]
    non_shared_Llsv_a_junc_coords = list(set(flattened_lsv_a_juncs).difference({shared_junc_coord_lsv_a}))
    LSV_A_minj = min(non_shared_Llsv_a_junc_coords)
    LSV_A_maxj = max(non_shared_Llsv_a_junc_coords)

    LSV_B_juncs = copy.copy(lsv_b['Junctions coords'])
    LSV_B_juncs = [LSV_B_juncs[x].split("-") for x in lsv_b["binary_indices"]]
    LSV_B_juncs = [map(int, x) for x in LSV_B_juncs]
    flattened_LSV_B_juncs = flatten_list(LSV_B_juncs)
    if len(set(flattened_LSV_B_juncs)) != 3:
        return "Neither"

    shared_junc_coord_LSV_B = list_duplicates(flattened_LSV_B_juncs)[0]
    non_shared_LSV_B_junc_coords = list(set(flattened_LSV_B_juncs).difference({shared_junc_coord_LSV_B}))
    LSV_B_minj = min(non_shared_LSV_B_junc_coords)
    LSV_B_maxj = max(non_shared_LSV_B_junc_coords)
    LSV_A_exons = set(lsv_a["Exons coords"])
    LSV_B_exons = set(lsv_b["Exons coords"])
    shared_exons_to_consider = list(LSV_A_exons.intersection(LSV_B_exons))
    shared_exons_to_consider = string_to_int_coords(shared_exons_to_consider)
    LSV_a_juncs = [[shared_junc_coord_lsv_a, LSV_A_minj], [shared_junc_coord_lsv_a, LSV_A_maxj]]
    LSV_b_juncs = [[shared_junc_coord_LSV_B, LSV_B_minj], [shared_junc_coord_LSV_B, LSV_B_maxj]]
    n_matched_exons_shared, matched_exons = count_shared_utilized_juncs(LSV_a_juncs, LSV_b_juncs,
                                                                        shared_exons_to_consider)

    LSV_A_binary_exons = get_binary_exons(lsv_a)
    LSV_B_binary_exons = get_binary_exons(lsv_b)
    if not isinstance(LSV_A_binary_exons[0], list) or not isinstance(LSV_B_binary_exons[0], list):
        return "Neither"  # alt 5' or 3'
    LSV_A_ref_exon = get_reference_exon(lsv_a)
    LSV_B_ref_exon = get_reference_exon(lsv_b)
    LSV_A_ref_exon_str = "-".join([str(coord) for coord in LSV_A_ref_exon])
    LSV_B_ref_exon_str = "-".join([str(coord) for coord in LSV_B_ref_exon])
    all_binary_exons = list()
    all_binary_exons.extend(LSV_A_binary_exons)
    all_binary_exons.extend(LSV_B_binary_exons)
    all_binary_exons_str = int_to_string_coords(all_binary_exons)
    all_binary_exons_str = set(all_binary_exons_str)
    LSV_A_binary_exons_str = ["-".join([str(coord) for coord in x]) for x in LSV_A_binary_exons]
    LSV_B_binary_exons_str = ["-".join([str(coord) for coord in x]) for x in LSV_B_binary_exons]
    if len(all_binary_exons_str) != 4:
        # not sure what these are
        if n_matched_exons_shared == 2:
            if str(matched_exons[0] != str(matched_exons[1])):  # make sure not same exon
                return "ME"
        else:
            return "Neither"
    # if LSV_A["Gene ID"] == "ENSMUSG00000042589" or LSV_B["Gene ID"] == "ENSMUSG00000042589":
    #     pdb.set_trace()
    # make sure there are the correct number of juncs:
    # 4 junctions, 2 pairs of which share reference LSV coord make 6 total unique
    # I think the ones that fail here are actually tandem cassettes...
    n_unique_coords = len(set(flattened_lsv_a_juncs).union(set(flattened_LSV_B_juncs)))
    if n_unique_coords != 6:
        if n_unique_coords == 4:
            alt_exons = all_binary_exons_str.difference(set(int_to_string_coords(matched_exons)))
            LSV_A_alt_exon = set(LSV_A_binary_exons_str).intersection(alt_exons)
            LSV_B_alt_exon = set(LSV_B_binary_exons_str).intersection(alt_exons)
            LSV_A_alt_exon = string_to_int_coords(list(LSV_A_alt_exon))[0]
            LSV_B_alt_exon = string_to_int_coords(list(LSV_B_alt_exon))[0]
            LSV_A_boundary = [max(LSV_A_alt_exon), max(LSV_A_ref_exon)]
            if min(LSV_A_boundary) <= max(LSV_B_alt_exon) <= max(LSV_A_boundary):
                return "Neither"  # almost looks like a ME!
            else:
                return "TC"
        else:
            return "Neither"
    else:
        return "Neither"


def count_shared_utilized_juncs(LSV_A_Juncs, LSV_B_Juncs, Exons):
    number_of_matched_exons_by_coord = 0
    matched_exons = list()
    for exon in Exons:
        if check_if_alt_exon_shared(LSV_A_Juncs[0], LSV_B_Juncs[0], exon):
            number_of_matched_exons_by_coord += 1
            matched_exons.append(exon)
        if check_if_alt_exon_shared(LSV_A_Juncs[1], LSV_B_Juncs[1], exon):
            number_of_matched_exons_by_coord += 1
            matched_exons.append(exon)
    return number_of_matched_exons_by_coord, matched_exons


def get_exons_containing_junc(Data, LSV_ID, Junc):
    """
    Given one end of a junction and an LSV ID, return
        all the exons in that LSV that the Junc falls into.

        retuns exon coords as list of ints [[#,#], [#,#], ...]
    """
    LSV = io_voila_caleb.get_lsv(Data, LSV_ID)
    exons_str = io_voila_caleb.get_exons(LSV)
    exons_int = string_to_int_coords(exons_str)
    matched_exons_int = io_voila_caleb.get_exons_containing_coord(exons_int, Junc)
    return matched_exons_int


def get_juncs_that_share_exon(data, lsv_id_sh, junc):
    """
    Given LSV_ID and a Junc coord, find which other juncs in the LSV
        land in the same exon.

        Junc cannot be within the LSV's reference exon.
    """
    matched_exons_int = io_voila_caleb.get_exons_containing_junc(data, lsv_id_sh, junc)
    matched_exons_str = int_to_string_coords(matched_exons_int)
    ref_exon_int = get_reference_exon(io_voila_caleb.get_lsv(data, lsv_id_sh))
    ref_exon_str = int_to_string_coords([ref_exon_int])[0]
    if ref_exon_str in matched_exons_str:
        raise RuntimeError("Oops, junc %s within %s" % (junc, lsv_id_sh))


def opposite_type(LSV_TYPE):
    """
    If "source", return "target".. vice versa
    """
    if LSV_TYPE == "source":
        return "target"
    elif LSV_TYPE == "target":
        return "source"
    else:
        pdb.set_trace()


def is_junc_in_lsv(lsv, junction_coord):
    """
    Is the junction in the LSV?
    :param lsv: lsv dict
    :param junction_coord: either "###-###" or [###,###]
    :return: Boolean
    """
    io_voila_caleb.check_is_lsv(lsv)
    if isinstance(junction_coord, list):
        if isinstance(junction_coord[0], int):
            junc_str = int_to_string_coords([junction_coord])[0]
    elif isinstance(junction_coord, str):
        if "-" not in junction_coord:
            raise ValueError("Expect '###-###', but instead got: %s" % junction_coord)
        else:
            junc_str = junction_coord
    else:
        raise RuntimeError("Unexpected Argument...")
    lsv_juncs = io_voila_caleb.get_juncs(lsv)
    if junc_str not in lsv_juncs:
        return False
    else:
        return True


def find_lsvs_with_junc(data, junc):
    """
    :param data: quick import structure
    :param junc: either "###-###" or [###,###]
    :return: dict {comp : [lsv, ids, containing junc] }
    """
    lsv_ids = []
    for comparison in data:
        lsvs = io_voila_caleb.get_lsvs_quickly(data, io_voila_caleb.get_LSV_IDs(data[comparison]), comparison)
        for lsv in lsvs:
            if is_junc_in_lsv(lsv, junc):
                lsv_ids.append(lsv["LSV ID"])
    return lsv_ids


def is_junc_connected_to_utilized_exon(Data,
                                       LSV_ID,
                                       Junction_Int,
                                       All_cas_like_things,
                                       Junc_maxPSI_dict,
                                       PSI_thresh=0.05):
    """
    Check if the provided junction is connected to a utilized exon.

        Basically, this function assumes that since you've provided the LSV ID,
        the LSV must be use the exon (I check to make sure of this, though...)

        Use the 'type' (source/target) opposite to the refence LSV's type to
        look for LSVs that include the exon at or above PSI_thresh.

        NOTE: if junction goes into another LSV, that LSV must be of the
            SAME type as the input LSV_ID to prove that it propegates the
            splicegraph in the same direction as LSV_ID. Otherwise, it
            could just be a cassette event (LSVs sharing an exclusion junction)

    """
    LSV = io_voila_caleb.get_lsv(Data, LSV_ID)
    if not is_junc_in_lsv(LSV, Junction_Int):
        raise RuntimeError("Please ensure junction is in LSV...")
    ref_chrm = io_voila_caleb.get_chr(LSV)
    ref_strand = io_voila_caleb.get_strand(LSV)
    Junction_str = int_to_string_coords([Junction_Int])[0]
    ref_junc_coord = ref_chrm + "_" + ref_strand + "_" + Junction_str
    if Junc_maxPSI_dict[ref_junc_coord] < PSI_thresh:
        print("Warning, junction %s isn't even used by LSV %s..." % (Junction_Int, LSV_ID))
        return False
    lsv_type_to_check = opposite_type(LSV["Reference_Type"])
    source_target_dict = get_sources_and_targets(Data, LSV_ID, Junction_Int)
    same_type_lsvs = source_target_dict[LSV["Reference_Type"]]
    for same_type in same_type_lsvs:
        LSV_to_check = io_voila_caleb.get_lsv(Data, same_type)
        ex_to_ch_str = io_voila_caleb.get_exons(LSV_to_check)
        ex_to_ch_str = set(list(ex_to_ch_str))
        ref_to_check_int = get_reference_exon(LSV_to_check)
        if is_junc_half_in_exon(Junction_Int, ref_to_check_int):
            if check_if_lsv_utilizes_intron(Data, same_type, 0.05):
                return True  # lets go with this is utilized...
            return "PROPOGATES"  # this means junction goes into LSV that propogates
            # the splicegraph same direction as LSV_ID... def utilized,
            # but could go on to earlier part of gene
    # now go on to look for LSVs that use exons your junction goes into..
    lsvs_to_check = source_target_dict[lsv_type_to_check]
    refLSV_of_interest_ex = get_reference_exon(LSV)
    if len(lsvs_to_check) == 0:
        return False
    redundant_lsv_found = False
    for lsv in lsvs_to_check:
        if check_if_lsv_utilizes_intron(Data, lsv, 0.05):
            return True  # lets go with this is utilized...
        LSV_to_check = io_voila_caleb.get_lsv(Data, lsv)
        # chr_strand_coords:max PSI
        ex_to_ch_str = io_voila_caleb.get_exons(LSV_to_check)
        ex_to_ch_str = set(list(ex_to_ch_str))  # remove dups
        # junction of interest may go into the lsv to check,
        # in which case this exon is def. 'utilized', but
        # very likely as a cassette... (share excl coords)
        # so, this function will just continue on looking
        # because this is redundant info (already know this
        # exon is utilized from the point of view of arg-provided junc)
        ref_to_check_int = get_reference_exon(LSV_to_check)
        if is_junc_half_in_exon(Junction_Int, ref_to_check_int):
            if lsv in All_cas_like_things:
                redundant_lsv_found = True
            # print "Could be cassette?", LSV_ID, lsv
            continue
        # else we need to see if the LSV utilizes the exon implicated by the
        # junction of interest.
        ref_to_check_str = int_to_string_coords([ref_to_check_int])[0]
        # already checked reference LSV
        ex_to_ch_str.remove(ref_to_check_str)
        ex_to_ch_int = string_to_int_coords(ex_to_ch_str)
        # get end of junc that isn't in refLSV of interest:
        # pdb.set_trace()
        non_ref_junc = get_non_ref_junc_coord(refLSV_of_interest_ex, Junction_Int)
        # get_non_ref_junc_coord returns a string if it can't figure out
        # which end of the junction is outside of the reference LSV's exon
        # for example 'NOVEL_INTRON_RETENTION'
        if isinstance(non_ref_junc, str):
            return False
        ex_to_ch_int = io_voila_caleb.get_exons_containing_coord(ex_to_ch_int, non_ref_junc)
        # make doubly sure that the junction of interest is hitting a non-ref
        # exon of some sort
        if len(ex_to_ch_int) < 1:
            # must be a funky case
            pdb.set_trace()
            return False
        if len(ex_to_ch_int) > 1:
            # overlapping?
            pdb.set_trace()
        # figure out which junction uses the exon
        juncs_to_ch = get_juncs_using_exon(LSV_to_check, ex_to_ch_int[0], Return_indices=True)
        lsv_juncs = io_voila_caleb.get_juncs(LSV_to_check)
        max_psi = 0
        for junc_ii in juncs_to_ch:
            this_junc = lsv_juncs[junc_ii]
            chrm = io_voila_caleb.get_chr(LSV_to_check)
            strand = io_voila_caleb.get_strand(LSV_to_check)
            junc_coord = chrm + "_" + strand + "_" + this_junc
            this_junc_max_psi = Junc_maxPSI_dict[junc_coord]
            max_psi = max(max_psi, this_junc_max_psi)
        if max_psi >= PSI_thresh:
            return True
        return False
    # if all LSVs cycled through, then non-utilized exon (or missing info in ttxt file)
    # or there was an LSV, but it is cassette-like (shares excl coord)
    if redundant_lsv_found:
        return True
    return False


def get_juncs_using_exon(LSV, Exon, Return_indices=True):
    """
    Given an Exon, return the junction(s) that go into it.

    Args:
        Return_indices: if True, return index of junction. If false,
            return the string '#-#'
    """
    io_voila_caleb.check_is_lsv(LSV)
    ref_exon = get_reference_exon(LSV)
    if are_exons_overlapping(ref_exon, Exon):
        pdb.set_trace()
    juncs_str = io_voila_caleb.get_juncs(LSV)
    juncs_int = string_to_int_coords(juncs_str)
    juncs_in_intron = list()
    junc_ii_in_intron = list()
    ii = 0
    for junc in juncs_int:
        Ju = get_non_ref_junc_coord(ref_exon, junc)
        if is_junc_in_exon(Ju, Exon):
            juncs_in_intron.append(junc)
            junc_ii_in_intron.append(ii)
        ii += 1
    if Return_indices:
        return junc_ii_in_intron
    return juncs_in_intron


def are_exons_overlapping(Exon1, Exon2):
    """
    If exons are overlapping, return True
    """
    if max(Exon1) > min(Exon2):
        if min(Exon1) < min(Exon2):
            return True
    if max(Exon2) > min(Exon1):
        if min(Exon2) < min(Exon1):
            return True
    return False


def get_sources_and_targets(Data, LSV_ID, Junction_int):
    """
    Given a Junction and LSV_ID, return the LSVs that have junctions
        pointing at the same (non-reference) exon.

        Returns dictionary{'sources': [lsvids], 'targets': [lsvids]}
    """
    res = dict()
    res["source"] = list()
    res["target"] = list()
    ref_lsv = io_voila_caleb.get_lsv(Data, LSV_ID)
    ref_ex_int = get_reference_exon(ref_lsv)
    non_ref_ju_int = get_non_ref_junc_coord(ref_ex_int, Junction_int)
    Gene_ExCoords_to_lsv_dict = Junc_to_LSVs(Data, LSV_ID, non_ref_ju_int)
    if len(Gene_ExCoords_to_lsv_dict) != 1:
        # if nove intron removal intron, skrew it!
        if check_if_novel_intron_lsv(Data, LSV_ID):
            return res
        pdb.set_trace()
    lsv_ids = Gene_ExCoords_to_lsv_dict.values()[0]
    lsv_ids.remove(LSV_ID)
    if len(lsv_ids) == 0:  # empty ..
        return res
    for lsv_id in lsv_ids:
        LSV = io_voila_caleb.get_lsv(Data, lsv_id)
        lsvtype = LSV["Reference_Type"]
        if lsvtype == "source":
            res["source"].append(lsv_id)
        else:
            res["target"].append(lsv_id)
    return res


def check_if_novel_intron_lsv(Data, LSV_ID, dPSI_Thresh=0.05):
    """
    Check if this LSV is one of those weird cases where there is a novel
        intron removal event within the ref exon itself.
        Often these are in the 3' UTR ...

        dPSI_Thresh= if threshold over which non-binary LSV junctions are
            tested for the above behavior.
    """
    io_voila_caleb.check_is_quick_import(Data)
    io_voila_caleb.check_is_lsv_id(LSV_ID)
    LSV = Data[Data.keys()[0]][LSV_ID]
    ref_ex = get_reference_exon(LSV)
    juncs = string_to_int_coords(LSV["Junctions coords"])
    if LSV.has_key("binary_indices"):
        bin_inds = LSV["binary_indices"]
        if is_junc_completely_within_exon(juncs[bin_inds[0]], ref_ex):
            return True
        if is_junc_completely_within_exon(juncs[bin_inds[1]], ref_ex):
            return True
    else:
        dpsis = io_voila_caleb.get_dpsis(LSV)
        if len(dpsis) != len(juncs):
            pdb.set_trace()
        abs_dpsis = [abs(x) for x in dpsis]
        for dpsi, junc in zip(abs_dpsis, juncs):
            if dpsi >= dPSI_Thresh:
                if is_junc_completely_within_exon(junc, ref_ex):
                    return True
            # if dpsi of this junction is unbalanced
            # for example, only 1 junction showing real changes
            # then the other junctions should be considered.
            if abs(sum(dpsis)) > 0.025:
                if is_junc_completely_within_exon(junc, ref_ex):
                    return True
    return False


def get_non_ref_junc_coord(Ref_Exon_Int, Junction_Int):
    """
    Given an exon, return the Junc coord that falls outside.
        returns "NOVEL_INTRON_REMOVAL" if both ends of junction
        are within the exon.
    """
    if is_junc_completely_within_exon(Junction_Int, Ref_Exon_Int):
        return "NOVEL_INTRON_REMOVAL"
    if min(Ref_Exon_Int) <= Junction_Int[0] <= max(Ref_Exon_Int):
        return Junction_Int[1]
    elif min(Ref_Exon_Int) <= Junction_Int[1] <= max(Ref_Exon_Int):
        return Junction_Int[0]
    else:
        pdb.set_trace()


def get_ref_junc_coord(Ref_Exon_Int, Junction_Int):
    """
    Given an exon, return the Junc coord that falls inside.
        returns "NOVEL_INTRON_REMOVAL" if both ends of junction
        are within the exon.
    """
    if is_junc_completely_within_exon(Junction_Int, Ref_Exon_Int):
        return "NOVEL_INTRON_REMOVAL"
    if min(Ref_Exon_Int) <= Junction_Int[0] <= max(Ref_Exon_Int):
        return Junction_Int[0]
    elif min(Ref_Exon_Int) <= Junction_Int[1] <= max(Ref_Exon_Int):
        return Junction_Int[1]
    else:
        pdb.set_trace()


def is_junc_in_exon(Junc, Exon):
    """
    Given a single junc int, return True
    if it is within the exon
    if min(Exon)<=Junc<=max(Exon): True
    """
    if min(Exon) <= Junc <= max(Exon):
        return True
    else:
        return False


def is_junc_half_in_exon(Junction_int, Exon):
    """
    Check if one coord of Junc's coords is between
     the exon coords.
    """
    count = 0
    if min(Exon) <= min(Junction_int) <= max(Exon):
        count += 1
    if min(Exon) <= max(Junction_int) <= max(Exon):
        count += 1
    if count != 1:
        return False
    return True


def is_junc_completely_within_exon(Junc, Exon):
    """
    Given: Junc : [#,#]
            Exon : [#,#]

            Return True if both ends of junction are
            contained within the exon. False otherwise.
    """
    if min(Exon) <= min(Junc) <= max(Exon):
        if min(Exon) <= max(Junc) <= max(Exon):
            return True
    return False


def Junc_to_LSVs(Data, LSV_ID, Junc):
    """
    Given LSV ID and Junction coord,
        return dictionary of GeneID_ExonCoords:[list of LSVs]
        whereby only those ExonCoords that contain the junc coordinate
        (not coordinate*s*) are in the dictionary.
    """
    if isinstance(Junc, list):
        raise RuntimeError("Received [] instead of # ...")
    io_voila_caleb.check_is_quick_import(Data)
    io_voila_caleb.check_is_lsv_id(LSV_ID)
    gene_id = LSV_ID.split(":")[0]
    exon_to_lsv_dict = GeneExons_to_LSVs(Data=Data, Gene_ID=gene_id)
    gene_exon_keys = exon_to_lsv_dict
    exons = list()
    for key in gene_exon_keys:
        exon_coords_str = key.split("_")[1]
        exon_coords_int = string_to_int_coords([exon_coords_str])[0]
        exons.append(exon_coords_int)
    matched_exons = io_voila_caleb.get_exons_containing_coord(exons, Junc)
    matched_exons_str = int_to_string_coords(matched_exons)
    matched_gene_exons_str = [gene_id + "_" + ex_str for ex_str in matched_exons_str]
    matched_exon_to_lsv_dict = {k: exon_to_lsv_dict.get(k, None) for k in matched_gene_exons_str}
    return matched_exon_to_lsv_dict


def GeneExons_to_LSVs(Data, Gene_ID):
    """
    Generate dictionary of GeneID_ExonCoords:[list of LSVs]
    """
    io_voila_caleb.check_is_quick_import(Data)
    comparisons = Data.keys()
    LSV_comparisons = io_voila_caleb.lookup_everywhere(Data,
                                                       Gene_ID,
                                                       save_lsv_structure_lookup=False,
                                                       print_bool=False)
    results = dict()
    for comparison in comparisons:
        LSVs = LSV_comparisons[comparison]
        for lsv_id in LSVs.keys():
            LSV = LSVs[lsv_id]
            exon_coords_str = copy.copy(LSV["Exons coords"])
            ref_exon_int = get_reference_exon(LSV)
            ref_exon_str = int_to_string_coords([ref_exon_int])[0]
            for exon_str in exon_coords_str:
                dict_key = Gene_ID + "_" + exon_str
                if dict_key not in results:
                    results[dict_key] = list()
                    results[dict_key].append(lsv_id)
                else:
                    results[dict_key].append(lsv_id)
                    # make sure lsv ids unique
                    results[dict_key] = list(set(results[dict_key]))
    return results


def get_exons_containing_coord(exons, coord):
    """
    Given list of exon coords and a single coord,
        return list of exon coords that contain the coord:
        min(exon)<=Coord<=max(exon)
    """
    exons_containing_coord = list()
    for exon in exons:
        if min(exon) <= coord <= max(exon):
            exons_containing_coord.append(exon)
    return exons_containing_coord


def get_non_reference_exons(LSV):
    ref_ex = get_reference_exon(LSV)
    ref_ex_str = int_to_string_coords([ref_ex])[0]
    all_ex = io_voila_caleb.get_exons(LSV)
    all_ex.remove(ref_ex_str)
    all_ex_int = string_to_int_coords(all_ex)
    return all_ex_int


def get_reference_exon(LSV):
    """
    Given LSV, return reference coord.
    """
    io_voila_caleb.check_is_lsv(LSV)
    exon = LSV["LSV ID"].split(":")[1]
    return string_to_int_coords([exon])[0]


def get_binary_exons(LSV):
    io_voila_caleb.check_is_lsv(LSV)
    if not LSV.has_key("binary_indices"):
        raise RuntimeError("LSV doesnt't look binary.")
    LSV_juncs = copy.copy(LSV['Junctions coords'])
    LSV_juncs = [LSV_juncs[x].split("-") for x in LSV["binary_indices"]]
    LSV_juncs = [map(int, x) for x in LSV_juncs]
    flattened_LSV_juncs = flatten_list(LSV_juncs)
    n_coord_ends = len(set(flattened_LSV_juncs))
    if n_coord_ends != 3:
        raise RuntimeError("LSV has %s insead of 3 coord ends" % n_coord_ends)
    shared_junc_coord_LSV = list_duplicates(flattened_LSV_juncs)[0]
    ref_exon = get_reference_exon(LSV)
    LSV_exons = string_to_int_coords(list(set(LSV["Exons coords"])))
    # if shared coord isn't in the reference exon
    if not (min(ref_exon) <= shared_junc_coord_LSV <= max(ref_exon)):
        # then it is alt 5' or 3'
        return ref_exon
    else:
        non_shared_LSV_junc_coords = list(set(flattened_LSV_juncs).difference({shared_junc_coord_LSV}))
        LSV_minj = min(non_shared_LSV_junc_coords)
        LSV_maxj = max(non_shared_LSV_junc_coords)
        binary_exons = list()
        binary_exons.extend(get_exons_containing_coord(LSV_exons, LSV_minj))
        binary_exons.extend(get_exons_containing_coord(LSV_exons, LSV_maxj))
        if len(set(int_to_string_coords(binary_exons))) == 1:
            # it is an alt 5' or 3'
            return list(set(int_to_string_coords(binary_exons)))[0]
        if len(set(int_to_string_coords(binary_exons))) != 2:
            pdb.set_trace()
            raise RuntimeError("Supposedly binary LSV has odd # (%s) of binary exons..." % len(binary_exons))
        return binary_exons


def string_to_int_coords(coords, coord_joiner="-"):
    """
    convert list of ['###-###'] to list of [[###,###]]
    :param coords:
    :return:
    """
    if not isinstance(coords, list):
        coords = [coords]
    if not isinstance(coords[0], str) or not coord_joiner in coords[0]:
        raise ValueError("Arguments are not formatted correctly.")
    coords = copy.copy(coords)
    coords = [x.split(coord_joiner) for x in coords]
    coords = [map(int, x) for x in coords]
    return coords


def int_to_string_coords(coords, coord_join="-"):
    """
    Given integer formatted coordinates, return as strings.
    :param coords: [[###,###], [###,###], ...]
    :param coord_join: what to join coords by?
    :return: ['###-###', '###-###', ...]
    """
    if not isinstance(coords, list):
        raise ValueError("Expected coords in [[###,###], [###,###], ...]")
    if not isinstance(coords[0], list):
        coords = [coords]
    if not isinstance(coords[0][0], int):
        raise ValueError("Should be ints in your coords here...")
    coord_strs = list()
    for coord in coords:
        str_coords = coord_join.join([str(co) for co in coord])
        coord_strs.append(str_coords)
    return coord_strs


def flatten_list(List):
    return [item for sl in List for item in sl]


def list_duplicates(seq):
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in seq if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def check_if_alt_exon_shared(alt1_junc, alt2_junc, ex):
    is_shared = 0
    if (min(ex)) <= alt1_junc[0] <= max(ex):
        is_shared += 1
    if (min(ex)) <= alt1_junc[1] <= max(ex):
        is_shared += 1
    if (min(ex)) <= alt2_junc[0] <= max(ex):
        is_shared += 1
    if (min(ex)) <= alt2_junc[1] <= max(ex):
        is_shared += 1
    if is_shared == 0:
        return False
    elif is_shared == 2:
        return True
    elif is_shared > 2:  # possible overlap of annotated exons
        # check if the junction coords define an shorter exon that makes sense
        check1 = abs(min(alt2_junc) - max(alt1_junc))
        check2 = abs(min(alt1_junc) - max(alt2_junc))
        if check1 > check2:
            junc_df_ex = [min(alt1_junc), max(alt2_junc)]
        else:
            junc_df_ex = [min(alt2_junc), max(alt1_junc)]
        if junc_df_ex == ex:
            print("stuck in loop ...")
            pdb.set_trace()
        if check_if_alt_exon_shared(alt1_junc, alt2_junc, junc_df_ex):
            print("recursively checking alt exon...")
            return True
        else:
            print("Oddity here..")
            return False
    else:
        # print "hmm ..."
        return False


def get_cassettes(Excl_dict):
    """
    Given output from gen_excl_dict(),
     return subsets of the dictionary that have:
        single LSV
        two (cassette)
        three plus
    """
    single = list()
    cassette = list()
    matched_excl_non_cassette = list()
    three_plus = list()
    n_complex = 0
    total = 0
    for chrom in Excl_dict.keys():
        for strand in Excl_dict[chrom].keys():
            for excl_junc_coords in Excl_dict[chrom][strand].keys():
                data = Excl_dict[chrom][strand][excl_junc_coords]
                data["excl_coords"] = excl_junc_coords
                n_lsvs = len(data["lsv_ids"])
                # if "ENSMUSG00000029761:34685963-34687173:target" in data["lsv_ids"]:
                #     pdb.set_trace()
                if n_lsvs == 1:
                    single.append(data)
                elif n_lsvs == 2:
                    alt1 = data['alt_coords'][0].split("-")
                    alt1 = [int(x) for x in alt1]
                    alt2 = data['alt_coords'][1].split("-")
                    alt2 = [int(x) for x in alt2]
                    if len(data['alt_exons']) != 1:
                        matched_excl_non_cassette.append(data)
                        total += n_lsvs
                        if len(data['alt_exons']) != 2:
                            pdb.set_trace()
                        continue
                    ex = data['alt_exons'][0]
                    if check_if_alt_exon_shared(alt1, alt2, ex):
                        cassette.append(data)
                    else:
                        pdb.set_trace()
                elif n_lsvs > 2:
                    n_complex += n_lsvs
                    three_plus.append(data)
                else:
                    print("Why are there 0 here?")
                    pdb.set_trace()
                total += n_lsvs
    print("Found %s single LSV, %s cassette-like, %s matched excl non-cassette, and %s (%s LSVs) multi-LSV events "
          "from %s total." % (len(single),
                              len(cassette),
                              len(matched_excl_non_cassette),
                              len(three_plus),
                              n_complex,
                              total))
    return single, cassette, matched_excl_non_cassette, three_plus


def gen_excl_dict(Binary_coord_data):
    """
    Given standard output from get_binary_coords(),
        generate a dictionary such that:

        chr -> strand -> lsv_id -> exclusion_junction_coord -> {
            "lsv_ids" = [list of LSV IDs],
            "alt_coords" = [list of #]
            TBD# "ref_coords" = ['#-#','#-#']
        }
    """
    excl_dict = dict()
    n_lsvs = 0
    for lsv_id in Binary_coord_data.keys():
        data = Binary_coord_data[lsv_id]
        chrom = data["chr"]
        strand = data["strand"]
        alt_junc_coord = data["alt_junc_coord"]
        alt_junc_coords = data["alt_junc_coords"]
        excl_junc_coord = data["excl_junc_coord"]
        excl_junc_coords = data["excl_junc_coords"]
        alt_exons = data["alt_exons"]
        if chrom not in excl_dict:
            excl_dict[chrom] = dict()
        if not excl_dict[chrom].has_key(strand):
            excl_dict[chrom][strand] = dict()
        if not excl_dict[chrom][strand].has_key(excl_junc_coords):
            excl_dict[chrom][strand][excl_junc_coords] = dict()
            excl_dict[chrom][strand][excl_junc_coords]["lsv_ids"] = list()
            excl_dict[chrom][strand][excl_junc_coords]["alt_coords"] = list()
            excl_dict[chrom][strand][excl_junc_coords]["alt_exons"] = list()
            # excl_dict[chrom][strand][alt_junc_coords]["ref_coords"] = list()
        excl_dict[chrom][strand][excl_junc_coords]["lsv_ids"].append(lsv_id)
        excl_dict[chrom][strand][excl_junc_coords]["alt_coords"].append(alt_junc_coords)
        excl_dict[chrom][strand][excl_junc_coords]["alt_exons"].extend(alt_exons)
        alts = excl_dict[chrom][strand][excl_junc_coords]["alt_exons"]
        unique_alts = [list(x) for x in set(tuple(x) for x in alts)]
        excl_dict[chrom][strand][excl_junc_coords]["alt_exons"] = unique_alts
        excl_dict[chrom][strand][excl_junc_coords]["chr"] = chrom
        excl_dict[chrom][strand][excl_junc_coords]["strand"] = strand
        n_lsvs += 1
        # excl_dict[chrom][strand][alt_junc_coords]["ref_coords"].append()
    print("From %s LSVs, extracted %s exclusion coordinates..." % (len(Binary_coord_data), n_lsvs))
    return excl_dict


def check_if_exons_overlap(Exon1, Exon2):
    if min(Exon2) < max(Exon1) < max(Exon2):
        return True
    if min(Exon1) < max(Exon2) < max(Exon1):
        return True
    return False


def get_binary_coords(Binary_Data,
                      Save_oddities=False,
                      Return_oddity_ids=False):
    """
    Given Binary LSV data, return dictionary of LSV IDs with the
        following keys

        chr: "Chr#"
        strand: "+/-"
        alt_junc_coord: #
        ref_coords : [start, end]
        excl_junc_coord: #
        excl_junc_coords: "#-#"
        alt_exons: [list of [#, #] possible exons]

        Return_oddity_ids: if True, also return LSV IDs that correspond to
            LSVs that, while binary-like, do not behave like a classic event.
                ex) both alt and excl junction coordinates starting from LSV reference,
                    and end with the same coordinate. These occur if the coord
                    in the reference LSV is different (different donor site).
        Save_oddities: if True, keep the oddities in the returned dict, otherwise,
            remove them,

        alt_exons is list of exons for whose coordinates the alt_junc_coord falls between

        Where, excl_junc_coord refers to the binary junction coordinate
            furthest from a given LSV's reference exon. Thus, if a given LSV
            is defining a cassette event, we expect the excl_junc_coord to
            match the given cassette's other reference exon LSV.
    """
    check_is_binary_lsv_data(Binary_Data)
    new_lsv_dict = dict()
    one_comparison = Binary_Data.keys()[0]
    lsv_dict = Binary_Data[one_comparison]
    lsv_ids = io_voila_caleb.get_LSV_IDs(lsv_dict)
    alt_ref_same_target = list()
    alt_exon_ov_ref = list()
    not_enough_exons = list()
    oddities = list()
    for lsv_id in lsv_ids:
        # if lsv_id == 'ENSMUSG00000033157:45742777-45742955:target':
        #     pdb.set_trace()
        lsv = lsv_dict[lsv_id]
        lsv_type = lsv["Reference_Type"]  # Source or Target?
        binary_indices = copy.copy(lsv["binary_indices"])
        exon_coords_str = copy.copy(lsv["Exons coords"])
        exon_coords_str = set(exon_coords_str)
        exon_coords = list()
        for exon in exon_coords_str:
            exon_coords.append([int(i) for i in exon.split("-")])
        jun_coords = copy.copy(lsv["Junctions coords"])
        jun_coords = [jun_coords[ii] for ii in binary_indices]
        ref_coords_str = lsv_id.split(":")[1].split("-")
        exon_coords_str.remove("-".join(ref_coords_str))
        # for exon in string_to_int_coords(exon_coords_str):
        #     if check_if_exons_overlap(string_to_int_coords(["-".join(ref_coords_str)])[0],exon):
        #         alt_exon_ov_ref.append(lsv_id)
        #         new_lsv_dict[lsv_id]["alt_exons"] = "ODDITY"
        #         continue
        ref_coords = [int(i) for i in ref_coords_str]
        jun_coords_int = [[int(c) for c in x.split("-")] for x in jun_coords]
        junc_coords1 = jun_coords_int[0]
        junc_coords2 = jun_coords_int[1]
        junc_coords1_diff = abs(junc_coords1[0] - junc_coords1[1])
        junc_coords2_diff = abs(junc_coords2[0] - junc_coords2[1])
        if junc_coords1_diff > junc_coords2_diff:  # if junc1 is wider
            alt_junc_coords_int = junc_coords2  # junc2 is inclusion/alt
            alt_junc_coords = jun_coords[1]
            excl_junc_coords_int = junc_coords1  # junc1 is exclusion
            excl_junc_coords = jun_coords[0]
        else:  # else junc2 is exclusion...
            alt_junc_coords_int = junc_coords1
            alt_junc_coords = jun_coords[0]
            excl_junc_coords_int = junc_coords2
            excl_junc_coords = jun_coords[1]
        try:
            new_lsv_dict[lsv_id] = dict()
            new_lsv_dict[lsv_id]["chr"] = lsv["chr"]
            new_lsv_dict[lsv_id]["strand"] = lsv["strand"]
            alt_junc_coord = alt_junc_coords_int[0]
            if alt_junc_coord >= min(ref_coords) and alt_junc_coord <= max(ref_coords):
                alt_junc_coord = alt_junc_coords_int[1]
            new_lsv_dict[lsv_id]["alt_junc_coord"] = alt_junc_coord
            new_lsv_dict[lsv_id]["alt_junc_coords"] = alt_junc_coords
            new_lsv_dict[lsv_id]["ref_coords"] = ref_coords
            excl_junc_coord = excl_junc_coords_int[0]
            if excl_junc_coord >= min(ref_coords) and excl_junc_coord <= max(ref_coords):
                excl_junc_coord = excl_junc_coords_int[1]
            new_lsv_dict[lsv_id]["excl_junc_coord"] = excl_junc_coord
            new_lsv_dict[lsv_id]["excl_junc_coords"] = excl_junc_coords
            if alt_junc_coord == excl_junc_coord:
                alt_ref_same_target.append(lsv_id)
            if len(exon_coords) < 2:
                not_enough_exons.append(lsv_id)
                new_lsv_dict[lsv_id]["alt_exons"] = "ODDITY"
                continue
            new_lsv_dict[lsv_id]["alt_exons"] = list()
            for exon in exon_coords:
                if alt_junc_coord >= min(exon) and alt_junc_coord <= max(exon):
                    new_lsv_dict[lsv_id]["alt_exons"].append(exon)
            if len(new_lsv_dict[lsv_id]["alt_exons"]) == 0:
                not_enough_exons.append(lsv_id)
                new_lsv_dict[lsv_id]["alt_exons"] = "ODDITY"
        except:
            print("Unknown Error...")
            pdb.set_trace()
    oddities.extend(alt_ref_same_target)
    oddities = list(set(oddities))
    if not Save_oddities:
        for lsvid in oddities:
            new_lsv_dict.pop(lsvid)
    if Return_oddity_ids:
        return new_lsv_dict, oddities
    print("%s non-classical events discarded from %s total, leaving %s" % (len(oddities),
                                                                           len(lsv_ids),
                                                                           len(new_lsv_dict)))
    return new_lsv_dict


def singly_unique(data,
                  threshold=0.05,
                  sig_dpsi_thresh=0.1,
                  sig_prob_thresh=0.95,
                  evaluate_introns=True,
                  comparisons=False,
                  stringent_prob=0,
                  must_reciprocate=True,
                  firstpdb=False,
                  secondpdb=False,
                  unblank_the_data=True):
    """
    Identify sigLSVs in the data that change in a pattern unique to a comparison.

    This algorithm is sensitive to the thresholds. Borderline cases will
        always be an issue. So, please check the results by eye to make sure
        the results pass the ole' "Yeah that checks out" intuition test.

    Arguments:
        data : Quick Import
        threshold (default 0.05) : dPSI w/in +/-[threshold] of 0 are converted to 0
        sig_dpsi_thresh (default 0.2) : dPSI thresh for ID'ing sig LSVs
        sig_prob_thresh (default 0.95) : Prob(dPSI) thresh for ID'ing sig LSVs
        evaluate_introns : True/False, should intronic LSVs be considered?
        comparisons (default False) : return list of comparisons corresponding to columns?
        stringent_prob (default 0) : convert dPSI to 0 if P(dPSI) < stringent_prob_thresh
            default 0 means no dPSI will be converted to 0 by this argument
            0.95 is very stringent, and will vastly reduce power to detect unique changes
                however, these unique changes will be highly supported by the RNA seq data!
        must_reciprocate (default True) : must the uniquely called LSV dPSIs sum to 0?
            (CAREFUL: BEHAVIOUR IF SET TO FALSE NOT TESTED YET)
            (not quite 'binary-like' b/c LSV with 4 junctions (2 pos, 2 neg) is reciprocating)
            after converting dPSI to 1,0,-1, checks that the sum for a uniquely called
                LSV is 0. For example, if it is a binary-like LSV, one junction goes up, then
                the other junction must go down. In conjunction with stringent_prob,
                this will very stringently identify uniquely changing LSVs of high confidence
                that make sense with respect to how splicing works on a biological level.
                Note: you will lose a binary LSV where one junction goes up, and the other
                    one doesn't go down enough or with enough probability to pass your threshes...
        firstpdb : debug? Boolean
        secondpdb : debug? Boolean
        unblank_the_data : Boolean
            If False, does not unblank the data, so you can play with the unblanked data if you want.

    Algorithm steps ::
     --(1)--
     --Convert dPSI to -1, 0, or 1
        if:
          dPSI < (-)Threshold, make it (-1)
          dPSI > (+)Threshold, make it ( 1)
          (-)Threshold <= dPSI <= (+)Threshold, make it 0
          prob(dPSI) < stringent_prob, make it 0

     --(2)--
     --Identify which junctions in which comparisons agree/disagree
        For each junction, check if all columns agree or not.
        - if all columns are the same, ignore the junction.
        - else, identify what the least common #(s) is(are), and then
        assign the column(s) name(s) to that row.

     --(3)--
     --Detect whether there is 1 disagreement that is singly unique to 1 comparison
        After going through all the junctions if:
         - at least one row is assigned a *single* column name
         - all rows agree on any *single* column name assignments
         - LSV reciprocates, if applicable (see must_reciprocate definition)
        Then that LSV is assigned to that column (it is unique to that comparison)

     --(4)-- TBD
     --Detect groups of comparisons that agree/disagree (TBD)
    if:
     - at least one row is assigned more than one column name
     - all rows agree on column name assignments
    Then that LSV is assigned to a combination of the comparisons
    eg: Comparison1_Comparison2 <- names will be sorted before joined

                |A vs B|A vs C  | all other comparisons ...
    Junction:   | dPSI | dPSI   | ...
            0   |   #  |  #     | ...
            1   |   #  |  #     | ...
            2   |   #  |  #     | ...
            ... |   ...|  ...   | ...

    Returns dictionary
        {SinglyUnique: {Comparisons: [LSV IDs ...]},
         cococombo_breaker: (Comparison_Combos: [LSV IDs ...])} <-TBD
    """
    non_red_sets, blanked_lsvs_dict, all_num_dpsis_data = non_redundant_set(data=data,
                                                                            cutoff_dpsi=0.1,
                                                                            cutoff_psi=1,
                                                                            save_blanked_structure=True,
                                                                            return_numdpsis_dat=True)
    print("Building non_red id dictionary ...")
    id_dict = non_redundant_id_dict(non_red_sets,
                                    cutoff_dpsi=0.1,
                                    cutoff_psi=1)
    print("Getting all Prob(dPSI) data ...")
    all_prob_data, comp_names = get_num_prob(data, True)
    sig_subset = io_voila_caleb.subset_significant(data,
                                                   sig_dpsi_thresh,
                                                   sig_prob_thresh,
                                                   keep_introns=evaluate_introns)
    all_uniques = list()
    singly_unique_sig = dict()
    non_sig = {"no_comp_assoc": []}
    cococombo_breaker = dict()
    opposite_dict = dict()
    not_changing = []
    non_recriprocating = []
    unsure = []
    sig_not_singly_so = []
    sig_subset_ids = dict()
    all_sig = []
    if comparisons:
        for c in comparisons:
            these_ids = io_voila_caleb.get_LSV_IDs(sig_subset[c])
            sig_subset_ids[c] = these_ids
            all_sig.extend(these_ids)
        comparisons_ii = [comp_names.index(x) for x in comparisons]
        all_sig = list(set(all_sig))
    else:
        for c in data.keys():
            sig_subset_ids[c] = io_voila_caleb.get_LSV_IDs(sig_subset[c])
        all_sig = io_voila_caleb.get_all_LSV_IDs(sig_subset)
        comparisons_ii = [comp_names.index(x) for x in data.keys()]
    thelsvs = all_num_dpsis_data.keys()
    indeces_at_10_percent = percent_through_list(thelsvs, 0.1)
    print("Assigning %s LSVs to groups ... " % (len(thelsvs)))
    i = 0.0
    for lsv_id in thelsvs:
        if i > 0.0 and i in indeces_at_10_percent:
            print(str(indeces_at_10_percent[i]) + "% of LSVs processed...")
        i += 1.0
        if firstpdb:
            if isinstance(firstpdb, str):
                if lsv_id == firstpdb:
                    pdb.set_trace()
            else:
                pdb.set_trace()
        num_dpsis_data = all_num_dpsis_data[lsv_id]
        num_probs_data = all_prob_data[lsv_id]
        if comparisons:
            num_dpsis_data = num_dpsis_data[:, comparisons_ii]
            num_probs_data = num_probs_data[:, comparisons_ii]
        # convert dPSIs to 1s and 0s
        num_probs_data[num_probs_data < stringent_prob] = 0
        num_probs_data[num_probs_data >= stringent_prob] = 1
        between_thresh_inclus = (num_dpsis_data >= -threshold) & (num_dpsis_data <= threshold)
        above_thresh = num_dpsis_data > threshold
        below_thresh = num_dpsis_data < -threshold
        num_dpsis_data[between_thresh_inclus] = 0
        num_dpsis_data[above_thresh] = 1
        num_dpsis_data[below_thresh] = -1
        # Account for stringent_prob_thresh
        num_dpsis_data = num_dpsis_data * num_probs_data
        # how many 1s, 0s, and -1s per junction?
        # Three arrays: -1s, 0s, and 1s
        # [j1, j2, ...] , [j1, j2, ...] , [j1, j2, ...]
        simplified_counts = counts_per_row(num_dpsis_data, [-1, 0, 1])
        # np.array(simplified_counts) = above col-concatenated, so first row is -1s,
        # second row is 0s, and last row is +1s
        # columns are junctions
        #
        # Which junctions have a n=1 for -1s, 0s, or 1s?
        sing_unique_calls = np.sum(np.array(simplified_counts) == 1, 0)
        if max(sing_unique_calls) == 1:
            if secondpdb:
                if isinstance(secondpdb, str):
                    if lsv_id == secondpdb:
                        pdb.set_trace()
                else:
                    pdb.set_trace()
            un_0s = np.where(np.array(simplified_counts) == 1)[0]  # row indices for 1s
            ov_0s = np.where(np.array(simplified_counts) == 1)[1]  # cols indices for 1s
            comp_iis = list()
            # get the comps corresponding to the uniquely called -1/1s
            for un_0, ov_0 in zip(un_0s, ov_0s):
                comp_i = np.where(num_dpsis_data[ov_0,] == [-1, 0, 1][un_0])[0][0]
                comp_iis.append(comp_i)
            # If only one comp implicated for the LSV
            if len(set(comp_iis)) == 1:
                unique_comp_i = list(set(comp_iis))[0]
                unique_comp = comp_names[comparisons_ii[unique_comp_i]]
                # if sig id for this comparison
                if lsv_id in sig_subset_ids[unique_comp]:
                    if must_reciprocate:
                        # sum junctions (cols) for each comparison
                        summed_juncs = np.sum(num_dpsis_data, axis=0)
                        if summed_juncs[unique_comp_i] != 0:
                            # Not a reciprocating LSV!! Skip it.
                            non_recriprocating.append(lsv_id)
                            continue
                    if unique_comp not in singly_unique_sig:
                        singly_unique_sig[unique_comp] = [lsv_id]
                        all_uniques.append(lsv_id)
                    else:
                        singly_unique_sig[unique_comp].append(lsv_id)
                        all_uniques.append(lsv_id)
                # not sig id for this comparison
                else:
                    # Maybe it is only sig for other comparison(s)
                    if lsv_id in all_sig:
                        sig_not_singly_so.append(lsv_id)
                    elif unique_comp not in non_sig:
                        non_sig[unique_comp] = [lsv_id]
                    else:
                        non_sig[unique_comp].append(lsv_id)
            else:
                # junctions seen changing at most in one comparison,
                # but more than one comparison implicated.
                # Thus, different junctions changing in different comparisons.
                # Could be non-reciprocating. If so, classify it as such.
                if False in (num_dpsis_data.sum(axis=0) == 0):
                    if lsv_id in all_sig:
                        non_recriprocating.append(lsv_id)
                    else:
                        non_sig["no_comp_assoc"].append(lsv_id)
                else:
                    if lsv_id in all_sig:
                        unsure.append(lsv_id)
                    else:
                        non_sig["no_comp_assoc"].append(lsv_id)

        # else if no -1s and no 1s
        elif simplified_counts[0].sum() == 0 and simplified_counts[2].sum() == 0:
            not_changing.append(lsv_id)
        else:
            # how many different junctions saw -1, 0, and 1?
            junc_counts = np.sum(simplified_counts, 1)
            neg1 = junc_counts[0]
            pos1 = junc_counts[2]
            # if multiple junctions agree on a LSVs comparisons:
            if neg1 == pos1:
                if neg1 == 0 and pos1 == 0:
                    print("unsure0")
                    pdb.set_trace()
                # could be an LSV shared by multiple comparisons..
                if must_reciprocate:
                    if False in (num_dpsis_data.sum(axis=0) == 0):
                        if lsv_id in all_sig:
                            non_recriprocating.append(lsv_id)
                        else:
                            non_sig["no_comp_assoc"].append(lsv_id)
                        continue
                comps_accrding_to_juncs = "immastring"
                dpsi_dir_agree = False
                opposites = False
                for row_ii in range(0, num_dpsis_data.shape[0]):
                    neg_comps = np.where(num_dpsis_data[row_ii,] < 0)
                    pos_comps = np.where(num_dpsis_data[row_ii,] > 0)
                    if neg_comps[0].shape[0] > 0 and pos_comps[0].shape[0] > 0:
                        # Comps that completely disagree on directionality of same juncs!
                        opposites = True
                    elif neg_comps[0].shape[0] == 0 and pos_comps[0].shape[0] == 0:
                        # then this junction never changes
                        continue
                    this_comps = np.where(abs(num_dpsis_data[row_ii,]) > 0)[0]
                    if not isinstance(comps_accrding_to_juncs, str):
                        dpsi_dir_agree = np.array_equal(this_comps, comps_accrding_to_juncs)
                        if not dpsi_dir_agree:
                            # Different comps use different junctions
                            break
                    elif this_comps.shape[0] != 0:
                        comps_accrding_to_juncs = this_comps
                # col_comps_ov_0 = np.where(num_dpsis_data > 0)[1]
                # col_comps_un_0 = np.where(num_dpsis_data < 0)[1]
                # # check if under and over 0 agree on which comparisons belong together
                # dpsi_dir_agree = np.array_equal(col_comps_ov_0, col_comps_un_0)
                if dpsi_dir_agree:
                    comp_iis = comps_accrding_to_juncs.tolist()
                    comps = [comp_names[comparisons_ii[comp_ii]] for comp_ii in comp_iis]
                    comps.sort()
                    comps = "_and_".join(comps)
                    if lsv_id in all_sig:
                        if opposites:
                            if comps in opposite_dict:
                                opposite_dict[comps].append(lsv_id)
                            else:
                                opposite_dict[comps] = [lsv_id]
                            continue
                        if comps in cococombo_breaker:
                            cococombo_breaker[comps].append(lsv_id)
                        else:
                            cococombo_breaker[comps] = [lsv_id]
                    elif comps in non_sig:
                        non_sig[comps].append(lsv_id)
                    else:
                        non_sig[comps] = [lsv_id]
                else:
                    if lsv_id in all_sig:
                        if False in (num_dpsis_data.sum(axis=0) == 0):
                            non_recriprocating.append(lsv_id)
                        else:
                            # Different comps use different junctions
                            unsure.append(lsv_id)
                    else:
                        non_sig["no_comp_assoc"].append(lsv_id)
            elif neg1 != pos1:
                # different number of -dPSI and +dPSI juncs cannot possibly reciprocate
                if lsv_id in all_sig:
                    non_recriprocating.append(lsv_id)
                else:
                    non_sig["no_comp_assoc"].append(lsv_id)
            else:
                print("unsure4")
                pdb.set_trace()

    print("Finished analyzing LSVs!")
    all_uniques = list(set(all_uniques))
    comp_names = list(singly_unique_sig.keys())
    comp_names.sort()
    nonrednetworks = dict()
    nonrednetworks["all_non_red_sets"] = non_red_sets
    nonrednetworks["unique_nonredsets"] = dict()
    nonrednetworks["cococombos"] = {x: [] for x in cococombo_breaker.keys()}
    nonrednetworks["opposites"] = {x: [] for x in opposite_dict.keys()}
    all_non_red_sets = non_red_sets['singles'].values()
    all_non_red_sets.extend(non_red_sets['twos'].values())
    all_non_red_sets.extend(non_red_sets['three_plus'].values())
    all_non_red_lsvs = non_red_sets['singles_all_lsvs']
    all_non_red_lsvs.extend(non_red_sets['twos_lsvs'])
    all_non_red_lsvs.extend(non_red_sets['three_plus_lsvs'])
    # n_sets = len(all_non_red_sets)
    n_lsvs = len(all_sig)
    summary_text = ""
    summary_text += "%s Sig LSVs using %s dPSI were grouped (or not)... ::\n" % (n_lsvs,
                                                                                 # n_sets,
                                                                                 sig_dpsi_thresh)
    if must_reciprocate:
        nonrednetworks["not_grouped"] = {"unsure": [],
                                         "sig_not_singly_so": [],
                                         "NonReciprocates": []}
    else:
        nonrednetworks["not_grouped"] = {"unsure": [],
                                         "sig_not_singly_so": []}
    # For each LSV identified as singly unique, generate a dict
    # that has the nonredsets that contain each LSV
    n_unique_lsvs = 0
    n_unique_sets = 0
    for comp in comp_names:
        nonrednetworks["unique_nonredsets"][comp] = []
        for sig_id in all_uniques:
            if sig_id in singly_unique_sig[comp]:
                nonrednetworks["unique_nonredsets"][comp].append(id_dict[sig_id])
        nonrednetworks["unique_nonredsets"][comp] = list(set(nonrednetworks["unique_nonredsets"][comp]))
        n_u = len(singly_unique_sig[comp])
        n_unique_lsvs += n_u
        n_nrn = len(nonrednetworks["unique_nonredsets"][comp])
        n_unique_sets += n_nrn
        summary_text += "%s: %s SinglyUnique LSVs, %s NonRedNetworks\n" % (comp, n_u, n_nrn)
    summary_text += "TOTAL Singly Uniques: %s LSVs, %s NonRedNetworks\n" % (n_unique_lsvs, n_unique_sets)
    n_combo_lsvs = 0
    n_combo_sets = 0
    for comp in cococombo_breaker.keys():
        n_ids = len(cococombo_breaker[comp])
        n_combo_lsvs += n_ids
        for sig_id in cococombo_breaker[comp]:
            if sig_id in all_sig:
                nonrednetworks["cococombos"][comp].append(id_dict[sig_id])
        nonrednetworks["cococombos"][comp] = list(set(nonrednetworks["cococombos"][comp]))
        n_nrn = len(nonrednetworks["cococombos"][comp])
        n_combo_sets += n_nrn
        summary_text += "===\n%s\n%s LSVs shared, %s NonRedNetworks\n" % (comp, n_ids, n_nrn)
    summary_text += "TOTAL combos: %s LSVs, %s NonRedNetworks\n" % (n_combo_lsvs, n_combo_sets)
    n_opp_lsvs = 0
    n_opp_sets = 0
    for comp in opposite_dict.keys():
        n_ids = len(opposite_dict[comp])
        n_opp_lsvs += n_ids
        for sig_id in opposite_dict[comp]:
            if sig_id in all_sig:
                nonrednetworks["opposites"][comp].append(id_dict[sig_id])
        nonrednetworks["opposites"][comp] = list(set(nonrednetworks["opposites"][comp]))
        n_nrn = len(nonrednetworks["opposites"][comp])
        n_opp_sets += n_nrn
        summary_text += "===\n%s\n%s LSVs show oppositivity, %s NonRedNetworks\n" % (comp, n_ids, n_nrn)
    summary_text += "TOTAL opposites: %s LSVs, %s NonRedNetworks\n" % (n_opp_lsvs, n_opp_sets)
    for sig_id in unsure:
        if sig_id in all_sig:
            nonrednetworks["not_grouped"]["unsure"].append(id_dict[sig_id])
    nonrednetworks["not_grouped"]["unsure"] = list(set(nonrednetworks["not_grouped"]["unsure"]))
    n_unsure_sets = len(nonrednetworks["not_grouped"]["unsure"])
    summary_text += "%s LSVs were not grouped because 'unsure,' %s NonRedNetworks\n" % (len(unsure), n_unsure_sets)
    for sig_id in sig_not_singly_so:
        if sig_id in all_sig:
            nonrednetworks["not_grouped"]["sig_not_singly_so"].append(id_dict[sig_id])
    nonrednetworks["not_grouped"]["sig_not_singly_so"] = list(set(nonrednetworks["not_grouped"]["sig_not_singly_so"]))
    n_sig_not_singly_so_sets = len(nonrednetworks["not_grouped"]["sig_not_singly_so"])
    summary_text += "%s LSVs were not grouped because 'sig_not_singly_so,' %s NonRedNetworks\n" % (
        len(sig_not_singly_so), n_sig_not_singly_so_sets)
    if must_reciprocate:
        for sig_id in non_recriprocating:
            if sig_id in all_sig:
                nonrednetworks["not_grouped"]["NonReciprocates"].append(id_dict[sig_id])
        nonrednetworks["not_grouped"]["NonReciprocates"] = list(set(nonrednetworks["not_grouped"]["NonReciprocates"]))
        n_lsvs_nonrec = len(non_recriprocating)
        n_nonrec_sets = len(nonrednetworks["not_grouped"]["NonReciprocates"])
        summary_text += "%s LSVs were not grouped because they're non-reciprocating, %s NonRedNetworks\n" % (
            n_lsvs_nonrec,
            n_nonrec_sets)
    print(summary_text)
    if unblank_the_data:
        io_voila_caleb.unimpute_lsv_data(data, blanked_lsvs_dict)
    results = {"singly_unique": singly_unique_sig,
               "cococombos": cococombo_breaker,
               "opposites": opposite_dict,
               "nonredsets": nonrednetworks,
               "summary": summary_text}
    return results


def counts_per_row(numpy_array: np.array,
                   elements_to_count):
    """


    :type numpy_array: numpy array duh
    :param numpy_array: np.array
    :param elements_to_count: list of elements to count, row-wise
    :return: list of arrays whereby each array corresponds to a row-wise count of how many
        times an element was seen. e.g.:

    numpy_array=
    array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.],
           [ 0.,  0.,  0., -1., -1., -1.,  0.,  0.]])

    elements_to_count=
            [-1, 0, 1]

    returns:
        [array([0, 0, 3]), array([8, 5, 5]), array([0, 3, 0])]
    """
    return [np.sum(numpy_array == x, axis=1) for x in elements_to_count]


def get_column_uniques(pandas_df, columns=False):
    """
    Given a binary pandas dataframe with only 1s and 0s,
        for each column, return a list of index names
        exclusive to that column. In other words, if only
        column Z has a 1 (all others cols have 0) in row X,
        then assign row X's index to column Z.

        Returns a dictionary {col_name: [list, of, indices]}
    """
    if columns:
        pandas_df_subset = pa.DataFrame(pandas_df, columns=columns)
    else:
        pandas_df_subset = pandas_df
    result_dict = dict()
    unique_ii_bools = pandas_df_subset.sum(axis=1) == 1
    pandas_df_unique_rows = pandas_df_subset[unique_ii_bools]
    for col_name in pandas_df_unique_rows.columns.tolist():
        are_ones_for_this_col = pandas_df_unique_rows[col_name] == 1
        unique_to_this_col = pandas_df_unique_rows[are_ones_for_this_col]
        result_dict[col_name] = unique_to_this_col.index.tolist()
    return result_dict


def sig_util_sets_to_binary_matrix(sig_util_non_red_sets):
    """
    Takes output from sig_utilized_non_red_sets and generates a
    numpy matrix of 0s and 1s. Each row is a set, each column is a comparison.

    Returns pandas dataframe with nonredsets as rownames and comparisons as col names.
    """
    all_SIGnonredsets = set()
    comparisons = sig_util_non_red_sets.keys()
    comparisons.sort()
    n_comp = len(comparisons)
    comparisons_i_dict = {comparison: comp_i for comparison, comp_i in zip(comparisons, range(0, n_comp))}
    for comparison in comparisons:
        this_comparisons_sets = sig_util_non_red_sets[comparison]
        all_SIGnonredsets = all_SIGnonredsets | this_comparisons_sets

    all_SIGnonredsets_list = list(all_SIGnonredsets)
    all_SIGnonredsets_list.sort()
    res_matrix = np.zeros((len(all_SIGnonredsets), len(comparisons)))
    for comparison in comparisons:
        this_comparisons_list = list(sig_util_non_red_sets[comparison])
        for nonred in this_comparisons_list:
            nonred_i = all_SIGnonredsets_list.index(nonred)
            if nonred_i == -1:  # then this set isn't significant in this comparison
                raise RuntimeError("Uhhh what? %s %s" % (nonred, comparison))
            res_matrix[nonred_i, comparisons_i_dict[comparison]] = 1
    dataframe = pa.DataFrame(res_matrix,
                             index=all_SIGnonredsets_list,
                             columns=comparisons)
    return dataframe


def sig_utilized_non_red_sets(Data,
                              Keep_introns=True,
                              CUTOFF_dPSI_nonred=0.2,
                              CUTOFF_PSI_nonred=1.0,
                              CUTOFF_dPSI_nonredset_hit=0.05,
                              CUTOFF_PROB_nonredset_hit=0,
                              dpsi_array_comparisons=False):
    """
    Generate a dict[comparison: sets] where the sets are the set of NonRedSets that
        the given comparison has at least one signififcant LSV ID in it.

        Also return non_red_dpsis[nonredset][lsvid]:dpsi_array

        CUTOFF_dPSI: what constitutes "significant"?
        CUTOFF_PROB: what constitutes "significant"?
        CUTOFF_PSI: what constitutes "significant"?
        dpsi_array_comparisons: list of comparisons to use for generating
            the dpsi_array results. Optional!
    """
    non_red_sets, blankdictids = non_redundant_set(data=Data,
                                                   cutoff_dpsi=CUTOFF_dPSI_nonred,
                                                   cutoff_psi=CUTOFF_PSI_nonred,
                                                   save_blanked_structure=True)
    id_dict = non_redundant_id_dict(non_red_sets,
                                    cutoff_dpsi=CUTOFF_dPSI_nonred,
                                    cutoff_psi=CUTOFF_PSI_nonred)
    # non_red_sets = set(id_dict.values())
    comparisons = Data.keys()
    results = {x: set() for x in comparisons}
    all_ids = io_voila_caleb.get_all_unique_lsv_ids(Data)
    n_a = len(all_ids)
    sig_dict_subset = io_voila_caleb.subset_significant(Data,
                                                        CUTOFF_dPSI_nonredset_hit,
                                                        CUTOFF_PROB_nonredset_hit,
                                                        Keep_introns)
    sig_ids = io_voila_caleb.get_all_unique_lsv_ids(sig_dict_subset)
    n_s = len(sig_ids)
    print("Identifying which NonRed sets are 'utilized' by sig LSVs in each comparison...")
    print("Out of %s total unique LSV IDs, %s are sig (dPSI>=%s,Prob>=%s)" % (n_a,
                                                                              n_s,
                                                                              CUTOFF_dPSI_nonredset_hit,
                                                                              CUTOFF_PROB_nonredset_hit))
    all_sets = get_all_nrsets(non_red_sets, Join="_")
    print("Building non-redundant set dPSI array dictionary...")
    nrsets_dpsis, dpsi_comps = non_redundant_dpsis(Data, all_sets, dpsi_array_comparisons)
    for comparison in comparisons:
        sig_ids = io_voila_caleb.get_LSV_IDs(sig_dict_subset[comparison])
        for sig_id in sig_ids:
            if sig_id in id_dict:
                results[comparison].add(id_dict[sig_id])
        # REVERT BLANKED DICT TO ORIGINAL STATE
        for lsv_id in blankdictids[comparison]:
            Data[comparison].pop(lsv_id)
    return results, nrsets_dpsis, dpsi_comps


def get_all_nrsets(NonRedSets, Join="_"):
    """
    Given value returned by Non_redundant_sets,
        return a list of all the non_redundant sets joined by underscore
    """
    loner_sets = NonRedSets["singles"].values()
    doublet_sets = NonRedSets["twos"].values()
    three_plus_sets = NonRedSets["three_plus"].values()
    sortListOLists(loner_sets)
    sortListOLists(doublet_sets)
    sortListOLists(three_plus_sets)
    all_sets = list()
    all_sets.extend(loner_sets)
    all_sets.extend(doublet_sets)
    all_sets.extend(three_plus_sets)
    all_sets_joined = list()
    for s in all_sets:
        s.sort()
        all_sets_joined.append(Join.join(s))
    return all_sets_joined


def sortListOLists(ListOLists, InPlace=True):
    if InPlace:
        for elem in ListOLists:
            elem.sort()
    else:
        copiedListOLists = copy.deepcopy(ListOLists)
        for elem in copiedListOLists:
            elem.sort()
        return copiedListOLists


def non_redundant_id_dict(data,
                          cutoff_dpsi=0.2,
                          cutoff_psi=1.0):
    """
    Given a quick import of data, generate a dict with:
        LSV_IDs -> partner_lsvs_sep_by_underscores
        e.g. 'ENSG00000004478:2906304-2906448:target' ->
            'ENSG00000004478:2906304-2906448:target_ENSG00000004478:2906304-2906448:source'

        Data: may be quick import of data OR value returned from Non_redundant_set().
        CUTOFF_dPSI: see Non_redundant_set() for details
        CUTOFF_PSI: see Non_redundant_set() for details
    """
    if check_is_non_red_set(data, True):
        non_red_sets = data
    else:
        non_red_sets = non_redundant_set(data=data,
                                         cutoff_dpsi=cutoff_dpsi,
                                         cutoff_psi=cutoff_psi)
    non_red_dict = dict()
    all_sets = get_all_nrsets(non_red_sets, Join="_")
    for nset in all_sets:
        nset_dict = {x: nset for x in nset.split("_")}
        non_red_dict.update(nset_dict)
    return non_red_dict


def non_redundant_dpsis(BlankedData,
                        NonRedSets,
                        Comparisons=False):
    """
    Given blanked quick import and NonRedSet,
        return:
        nrset_numdpsis_dict[nonredsetid] : lsvid : numdPSIs array

        Comparisons:
            Optional! If list of names provided, subset the
            BlankedData to only include those comparisons.
    """
    if isinstance(Comparisons, str):
        comparisons = list(Comparisons)
    elif isinstance(Comparisons, list):
        comparisons = Comparisons
    else:
        comparisons = BlankedData.keys()
        comparisons.sort()
    BlankedData = {k: BlankedData[k] for k in comparisons}
    all_numdpsis, comparisons = get_num_d_psi(BlankedData, True)
    nrset_numdpsis = dict()
    for nrset in NonRedSets:
        nrset_numdpsis[nrset] = dict()
        lsv_ids = nrset.split("_")
        for lsv_id in lsv_ids:
            nrset_numdpsis[nrset][lsv_id] = all_numdpsis[lsv_id]
    return nrset_numdpsis, comparisons


def find_set_partners(connected_sets, lsv_id, sig_ids=False):
    """
    Given data returned from Non_redundant_set and a query LSV_ID,
        return the LSV's partenr(s) that share a utilized junction.

        Arguments:
            connected_sets: return value from Non_redundant_set
            lsv_id: LSV_ID: query LSV to get partnrs for
            sig_ids: OPTIONAL list/set
                If provided, this funciton will only return partner LSV(s)
                that are in this list/set. For example, if you only
                want to identify "significant" LSV partners...
    """
    if sig_ids:
        sig_ids = list(sig_ids)
    if lsv_id in connected_sets['singles_all_lsvs']:
        return "isolated"
    elif lsv_id in connected_sets['twos_lsvs']:
        pairs = connected_sets['twos'].values()
        for pair in pairs:
            if lsv_id in pair:
                partner = copy.copy(pair)
                partner.remove(lsv_id)
                if sig_ids:
                    if partner[0] in sig_ids:
                        return partner
                    else:
                        return "no_sig_partner"
                else:
                    return partner
    elif lsv_id in connected_sets['three_plus_lsvs']:
        complexes = connected_sets['three_plus'].values()
        for comp in complexes:
            if lsv_id in comp:
                partners = copy.copy(comp)
                partners.remove(lsv_id)
                if sig_ids:
                    sig_partners = list(set(sig_ids).intersection(set(partners)))
                    if len(sig_partners) > 0:
                        return partners
                    else:
                        return "no_sig_partners"
                else:
                    return partners
    else:
        print(lsv_id, "wasn't found in sets...")


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
    blanked_lsvs_dict = io_voila_caleb.impute_missing_lsvs(data=data, in_place=True, warnings=False)
    print("Finished filling in the gaps, running non-redundant algorithm...")
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
    print("Accounting for LSVs from different genes that overlap according to genomic coordinates ...")
    by_type = most_lsvs_same_gene(Connected_LSVs=nr_connected_lsvs)
    if not save_blanked_structure:
        print("Reverting Data to original state....")
        for comparison in data.keys():
            for lsv_id in blanked_lsvs_dict[comparison]:
                data[comparison].pop(lsv_id)
    print("Done!!!\n\n\n")
    s = len(by_type["singles"])
    ns = len(by_type["singles_all_lsvs"])
    d = len(by_type["twos"])
    nd = len(by_type["twos_lsvs"])
    tp = len(by_type["three_plus"])
    ntp = len(by_type["three_plus_lsvs"])
    total = s + d + tp
    total_lsvs = ns + nd + ntp
    print(total_lsvs, "total LSVs met junction utilization cutoffs...")
    print(total, "Total non-redundant splicing sets")
    print("%s singles (%s LSVs), %s doubles (%s LSVs), %s three_pluses (%s LSVs) ..." % (s, ns, d, nd, tp, ntp))
    if return_numdpsis_dat:
        if save_blanked_structure:
            return by_type, blanked_lsvs_dict, nr_numdpsis
        return by_type, nr_numdpsis
    else:
        if save_blanked_structure:
            return by_type, blanked_lsvs_dict
        return by_type


def most_lsvs_same_gene(Connected_LSVs,
                        Return_shared_junc_genes=False):
    """
    Some lsvs that are connected by a junction turn out to
        to be due to overlapping genes. Thus, this function teases
        out how which connected LSVs correspond to connections in
        the same gene. Removes LSVs that are redundant.

        Return_shared_junc_genes: if True, return list of junctions that
            multiple LSVs (in more than 1 gene) share in the Connected
            LSV network.
    """

    lengths = [len(Connected_LSVs[x]) for x in Connected_LSVs.keys()]
    lengths = np.array(lengths)
    singles = np.where(lengths == 1)[0]
    singles = singles.tolist()
    singles_lsvs = list()
    for single in singles:
        singles_lsvs.extend(Connected_LSVs[Connected_LSVs.keys()[single]])
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
        connected = Connected_LSVs[Connected_LSVs.keys()[two_plus]]
        gene_ids = [x.split(":")[0] for x in connected]
        if len(set(gene_ids)) > 1:
            genes_share_junc.append(Connected_LSVs.keys()[two_plus])
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
    results["singles"] = {Connected_LSVs.keys()[k]: Connected_LSVs[Connected_LSVs.keys()[k]] for k in singles}
    results["singles_all_lsvs"] = singles_lsvs
    results["twos"] = {Connected_LSVs.keys()[k]: Connected_LSVs[Connected_LSVs.keys()[k]] for k in twos}
    results["twos_lsvs"] = two_lsvs
    results["three_plus"] = {Connected_LSVs.keys()[k]: Connected_LSVs[Connected_LSVs.keys()[k]] for k in three_plus}
    results["three_plus_lsvs"] = three_lsvs
    if Return_shared_junc_genes:
        dict_res = dict()
        dict_res["types"] = results
        dict_res["genes_share_junc"] = genes_share_junc
        return dict_res
    return results


def most_common(L):
    """
    return most common element of list

    http://stackoverflow.com/a/1520716
    """
    # get an iterable of (item, iterable) pairs
    SL = sorted((x, i) for i, x in enumerate(L))
    # print 'SL:', SL
    groups = itertools.groupby(SL, key=operator.itemgetter(0))

    # auxiliary function to get "quality" for an item
    def _auxfun(g):
        item, iterable = g
        count = 0
        min_index = len(L)
        for _, where in iterable:
            count += 1
            min_index = min(min_index, where)
        # print 'item %r, count %r, minind %r' % (item, count, min_index)
        return count, -min_index

    # pick the highest-count/earliest item
    return max(groups, key=_auxfun)[0]


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
    all_lsvs = io_voila_caleb.get_shared_slv_ids(data)
    junc_lsv_dicts = get_junc_lsv_dicts(data,
                                        cutoff_dpsi=Cutoff_dPSI,
                                        cutoff_psi=Cutoff_PSI,
                                        return_numdpsis=ret_numpdsis_data)
    junc_to_lsv = junc_lsv_dicts["junc_lsv_dict"]
    lsv_to_junc = junc_lsv_dicts["lsv_junc_dict"]
    sig_juncs = junc_to_lsv.keys()
    print("Recursively processing %s sig junctions to build connectivity graph ..." % (len(sig_juncs)))
    indeces_at_10_percent = percent_through_list(sig_juncs, 0.01)
    i = 0.0
    n_to_start = len(sig_juncs)
    while len(sig_juncs) > 0:
        if i > 0.0 and i in indeces_at_10_percent:
            print(str(indeces_at_10_percent[i]) + "% of sig juncs processed...")
        next_junc = sig_juncs.pop()
        conn_lsvs = connected_lsvs(all_lsvs, sig_juncs, junc_to_lsv, lsv_to_junc, next_junc)
        master_junc_dict[next_junc] = conn_lsvs
        i = n_to_start - len(sig_juncs)
    if ret_numpdsis_data:
        return master_junc_dict, junc_lsv_dicts["numdPSIs_data"]
    return master_junc_dict


def connected_lsvs(All_lsvs, Sig_juncs, Junc_to_LSV, LSV_to_Junc, Start_junc):
    """
    Return all the LSVs connected to the start junc, recursively.
    """
    conn_lsvs = Junc_to_LSV[Start_junc]
    associated_juncs = list()
    associated_lsvs = list()
    # if Start_junc == 'chr11_+_8705628-8706265':
    #     pdb.set_trace()
    # if Start_junc == 'chr11_+_8705628-8706370':
    #     pdb.set_trace()
    for found_lsv in conn_lsvs:
        if found_lsv not in All_lsvs:
            continue
        else:
            associated_lsvs.append(found_lsv)
            All_lsvs.remove(found_lsv)
            found_lsv_juncs = LSV_to_Junc[found_lsv]
            found_lsv_juncs = found_lsv_juncs.tolist()
            associated_juncs.extend(found_lsv_juncs)
    associated_juncs = [x for x in associated_juncs if x in Sig_juncs]  # only keep unseen juncs
    for found_junc in associated_juncs:
        if found_junc not in Sig_juncs:
            continue
        else:
            Sig_juncs.remove(found_junc)
            recursive_call = connected_lsvs(All_lsvs,
                                            Sig_juncs,
                                            Junc_to_LSV,
                                            LSV_to_Junc,
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
    if io_voila_caleb.check_is_lsv_dict(Data, da_bool=True):
        master_junc_weight_dict = get_junc_weights(Data)
    elif io_voila_caleb.check_is_quick_import(Data):
        shared_ids = io_voila_caleb.get_shared_slv_ids(Data)
        comp_to_junc_weights = dict()
        all_comparisons = Data.keys()
        for comparison in all_comparisons:
            share_dict = io_voila_caleb.lsv_dict_subset(Data[comparison], shared_ids)
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
    io_voila_caleb.check_is_quick_import(Data)
    comparisons = Data.keys()
    shared_ids = io_voila_caleb.get_shared_slv_ids(Data)
    shared_dicts = dict()
    for comparison in comparisons:
        share_dict = io_voila_caleb.lsv_dict_subset(Data[comparison], shared_ids)
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
    io_voila_caleb.check_is_lsv_dict(Data)
    lsv_ids = io_voila_caleb.get_LSV_IDs(Data)
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
    io_voila_caleb.check_is_quick_import(data)
    print("Getting all dPSI data ...")
    numdPSIs_data = get_num_d_psi(data)
    print("Getting all PSI data ...")
    numPSIs_data = get_num_psi(data)
    print("Getting all junction coordinates ...")
    lsv_junc_dict = get_junc_coords(data)
    print("Running find_binary_LSV_IDs() ...")
    sig_junc_dict = find_binary_lsv_ids(num_d_psi=numdPSIs_data,
                                        this_num_psi=numPSIs_data,
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
        res["numdPSIs_data"] = numdPSIs_data
    return res


def get_junc_coords(Data):
    """
    Given LSV dict info (quick import, or LSV dict),
     return dict of lsv_ids pointing at array of junction loci:
        'chr_strand_start_end'
    """
    if io_voila_caleb.check_is_quick_import(Data, the_bool=True):
        junc_dict = dict()
        for comp in Data:
            junc_dict.update(get_junc_coords(Data[comp]))
        return junc_dict
    io_voila_caleb.check_is_lsv_dict(Data)
    junc_dict = dict()
    lsv_ids = io_voila_caleb.get_LSV_IDs(Data)
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


def get_num_d_psi(data,
                  return_comparisons=False,
                  use_binary_index_info=None):
    """
    Given dictionary of LSV dictionaries return numpy arrays
        of all dPSIs from all comparisons for each LSV. Such that
        each LSV ID in the dictionary is pointing at an array that
        looks like this:

                |A vs B|A vs C  | all other comparisons ...
    Junction:   | dPSI | dPSI   | ...
            0   |   #  |  #     | ...
            1   |   #  |  #     | ...
            2   |   #  |  #     | ...
            ... |   ...|  ...   | ...

    Arguments:
        data: quick imp
        return_comparisons: Boolean. If True, also return a list
            of the comparison names, in the same order as the columns
            of the numpy array.
        use_binary_index_info: None, "closer", or "further"
            If the LSV dictionaries were returned by get_binary_LSVs(),
            then they will have a "binary_indices" key that points at
            the binary indices. Use this info to return only the dPSI
            value that corresponds to the closer or further junction
            from the reference exon.

    NOTES:
        - only those LSV IDs that are shared by all LSV dictionaries
        in Data are evaluated. LSV IDs that are unique to or missing from
        any of the LSV Dictionaries are not returned by this function.
        - columns are sorted by name


    Returns the numpy array and a list of LSV_IDs
    """
    binary_index = "cow"  # stupid pep
    if use_binary_index_info:
        if not use_binary_index_info == "closer" and not use_binary_index_info == "further":
            raise ValueError("Use_binary_index_info needs to be 'closer' or 'further' if provided")
    if use_binary_index_info == "closer":
        binary_index = 0
    if use_binary_index_info == "further":
        binary_index = 1
    io_voila_caleb.check_is_quick_import(data)
    comparison_dict = dict()
    for comparison in data.keys():
        d_psis = listdPSI(data[comparison])
        comparison_dict[comparison] = d_psis
    union_of_lsv_ids = io_voila_caleb.get_shared_slv_ids(data)
    lsv_dict = dict()
    comparisons = list(comparison_dict.keys())
    comparisons.sort()  # alphabatized
    for lsv_id in list(union_of_lsv_ids):
        list_of_dPSIs = list()
        for comparison in comparisons:
            d_psis = comparison_dict[comparison][lsv_id]
            if use_binary_index_info:
                binary_i = data[comparison][lsv_id]["binary_indices"][binary_index]
                d_psis = d_psis[binary_i]
            list_of_dPSIs.append(d_psis)
        lsv_dict[lsv_id] = np.array(list_of_dPSIs).T
    if return_comparisons:
        return lsv_dict, comparisons
    return lsv_dict


def get_num_prob(data,
                 return_comparisons=False,
                 use_binary_index_info=False):
    """
    Given dictionary of LSV dictionaries return numpy arrays
        of all P(dPSIs) from all comparisons for each LSV. Such that
        each LSV ID in the dictionary is pointing at an array that
        looks like this:

                |A vs B|A vs C  | all other comparisons ...
    Junction:   | P(dPSI) | P(dPSI)   | ...
            0   |   #     |     #     | ...
            1   |   #     |     #     | ...
            2   |   #     |     #     | ...
            ... |   ...   |     ...   | ...
    Arguments:
        data: Quick Import structure
        return_comparisons: Boolean. If True, also return a list
            of the comparison names, in the same order as the columns
            of the numpy array.
        use_binary_index_info: None, "closer", or "further"
            If the LSV dictionaries were returned by get_binary_LSVs(),
            then they will have a "binary_indices" key that points at
            the binary indices. Use this info to return only the P(dPSI)
            value that corresponds to the closer or further junction
            from the reference exon.

    NOTES:
        - only those LSV IDs that are shared by all LSV dictionaries
        in Data are evaluated. LSV IDs that are unique to or missing from
        any of the LSV Dictionaries are not returned by this function.
        - columns are sorted by name


    Returns the numpy array and a list of LSV_IDs
    """
    # stupid pep rules
    binary_index = "cow"
    if use_binary_index_info:
        if not use_binary_index_info == "closer" and not use_binary_index_info == "further":
            raise ValueError("Use_binary_index_info needs to be 'closer' or 'further' if provided")
        if use_binary_index_info == "closer":
            binary_index = 0
        elif use_binary_index_info == "further":
            binary_index = 1
        else:
            raise RuntimeError("I don't know what to do here")
    io_voila_caleb.check_is_quick_import(data)
    comparison_dict = dict()
    for nup_comparison in data.keys():
        nu_probs = list_probs(data[nup_comparison])
        comparison_dict[nup_comparison] = nu_probs
    union_of_lsv_ids = io_voila_caleb.get_shared_slv_ids(data)
    lsv_dict = dict()
    comparisons = list(comparison_dict.keys())
    comparisons.sort()  # alphabatized
    for lsv_id in list(union_of_lsv_ids):
        list_of_probs = list()
        for nup_comparison in comparisons:
            probs = comparison_dict[nup_comparison][lsv_id]
            if use_binary_index_info:
                binary_i = data[nup_comparison][lsv_id]["binary_indices"][binary_index]
                probs = probs[binary_i]
            list_of_probs.append(probs)
        lsv_dict[lsv_id] = np.array(list_of_probs).T
    if return_comparisons:
        return lsv_dict, comparisons
    return lsv_dict


def list_probs(lsv_dict):
    """
    Return dictionary of LSV IDs pointing at list of P(dPSIs):
    (this col
    not in dict) list:
    Junction:   | Prob |
            0   |   #  |
            1   |   #  |
            2   |   #  |
            ... |   ...|
    """
    io_voila_caleb.check_is_lsv_dict(lsv_dict)
    lsv_to_prob_dict = dict()
    lsvs = io_voila_caleb.get_LSV_IDs(lsv_dict)
    # Extract dPSI from each LSV, using cutoff.
    for lsv_id in lsvs:
        lsv_to_prob_dict[lsv_id] = io_voila_caleb.get_probs(lsv_dict[lsv_id])
    return lsv_to_prob_dict


def listdPSI(LSV_dict):
    """
    Return dictionary of LSV IDs pointing at list of dPSIs:
    (this col
    not in dict) list:
    Junction:   | dPSI |
            0   |   #  |
            1   |   #  |
            2   |   #  |
            ... |   ...|
    """
    io_voila_caleb.check_is_lsv_dict(LSV_dict)
    lsv_to_dPSI_dict = dict()
    lsvs = io_voila_caleb.get_LSV_IDs(LSV_dict)
    # Extract dPSI from each LSV, using cutoff.
    for lsv_id in lsvs:
        lsv_to_dPSI_dict[lsv_id] = io_voila_caleb.get_dpsis(LSV_dict[lsv_id])
    return lsv_to_dPSI_dict


def get_num_psi(data,
                return_comparisons=False,
                use_binary_index_info=None):
    """
    Given dictionary of LSV dictionaries return numpy array
        of all PSIs from all conditions for each LSV. Such that
        each LSV ID in the dictionary is pointing at an array that
        looks like this:

                |   A  |  B     | all other conditions ...
    Junction:   |  PSI |  PSI   | ...
            0   |   #  |  #     | ...
            1   |   #  |  #     | ...
            2   |   #  |  #     | ...
            ... |   ...|  ...   | ...

    Arguments:
        return_comparisons: if True, also return a list of conditions
            found in the same order as the numpy array columns.
        use_binary_index_info: None, "closer", or "further"
            If the LSV dictionaries were returned by get_binary_LSVs(),
            then they will have a "binary_indices" key that points at
            the binary indices. Use this info to return only the PSI
            value that corresponds to the closer or further junction
            from the reference exon.

    NOTE: only those LSV IDs that are shared by all LSV dictionaries
        in Data are evaluated. LSV IDs that are unique to or missing from
        any of the LSV Dictionaries are not returned by this function.

    Returns the numpy array and a list of LSV_IDs
    """
    if use_binary_index_info:
        if not use_binary_index_info == "closer" and not use_binary_index_info == "further":
            raise ValueError("Use_binary_index_info needs to be 'closer' or 'further' if provided")
    if use_binary_index_info == "closer":
        binary_index = 0
    if use_binary_index_info == "further":
        binary_index = 1
    comparison = "thing to make PyCharm happy"
    binary_index = "thing to make PyCharm happy"
    io_voila_caleb.check_is_quick_import(data)
    condition_dict = dict()
    lsv_id_lists = list()
    for comparison in data.keys():
        lsv_dict = data[comparison]
        cond_1_name = lsv_dict["condition_1_name"]
        cond_2_name = lsv_dict["condition_2_name"]
        PSIs = list_psi(lsv_dict)
        cond_1_PSIs = PSIs[cond_1_name]
        cond_2_PSIs = PSIs[cond_2_name]
        condition_dict[cond_1_name] = cond_1_PSIs
        condition_dict[cond_2_name] = cond_2_PSIs
    union_of_lsv_ids = io_voila_caleb.get_shared_slv_ids(data)
    lsv_dict = dict()
    conditions = list(condition_dict.keys())
    conditions.sort()
    for lsv_id in list(union_of_lsv_ids):
        # if lsv_id == "ENSG00000173744:228395807-228395926:target":
        #     pdb.set_trace()
        list_of_PSIs = list()
        for condition in conditions:
            PSI = condition_dict[condition][lsv_id]
            if use_binary_index_info:
                binary_i = data[comparison][lsv_id]["binary_indices"][binary_index]
                PSI = PSI[binary_i]
            list_of_PSIs.append(PSI)
        lsv_dict[lsv_id] = np.array(list_of_PSIs).T
    if return_comparisons:
        return lsv_dict, conditions
    return lsv_dict


def list_psi(lsv_dict):
    """
    Return dictionary of conditions pointing at dictionary of
        LSV IDs pointing at list of PSIs:

    condition_* points at:
    (this col
    not in dict) list:
    Junction:   | PSI  |
            0   |   #  |
            1   |   #  |
            2   |   #  |
            ... |   ...|
    """
    io_voila_caleb.check_is_lsv_dict(lsv_dict)
    cond_1_name = lsv_dict["condition_1_name"]
    cond_2_name = lsv_dict["condition_2_name"]
    condtion_dict = dict()
    condtion_dict[cond_1_name] = dict()
    condtion_dict[cond_2_name] = dict()
    lsvs = io_voila_caleb.get_LSV_IDs(lsv_dict)
    # Extract dPSI from each LSV, using cutoff.
    for lsv_id in lsvs:
        PSIs = io_voila_caleb.get_psis(lsv_dict[lsv_id])
        cond_1_psi = PSIs[0]
        cond_2_psi = PSIs[1]
        condtion_dict[cond_1_name][lsv_id] = cond_1_psi
        condtion_dict[cond_2_name][lsv_id] = cond_2_psi
    return condtion_dict


def check_is_non_red_set(NonRedSet, Bool=False):
    """
    Check that NonRedSet is, in fact, value returned by non_red_set()
    """
    if not isinstance(NonRedSet, dict):
        if not Bool:
            raise ValueError("Expected a NonRedSet, instead object is of type: "
                             + str(type(NonRedSet)))
        else:
            return False
    for key in ['twos_lsvs',
                'singles',
                'three_plus',
                'three_plus_lsvs',
                'singles_all_lsvs',
                'twos']:
        if key not in NonRedSet:
            if not Bool:
                raise ValueError("Expected a NonRedSet, but missing key: %s" % (key))
            else:
                return False
    return True


def check_is_binary_lsv_data(data,
                             thisbool=False):
    """
    """
    io_voila_caleb.check_is_quick_import(data)
    io_voila_caleb.check_lsv_ids_all_shared(data)
    lsv_dict = data[list(data.keys())[0]]
    random_lsv_id = io_voila_caleb.get_LSV_IDs(lsv_dict)[0]
    if "binary_indices" in lsv_dict[random_lsv_id]:
        if thisbool:
            return True
    elif thisbool:
        return False
    else:
        raise ValueError("Provided LSV dictionary data didn't have binary_indices")


def get_junc_index(data, lsv, junc_coord):
    io_voila_caleb.check_is_quick_import(data)
    io_voila_caleb.check_is_lsv(lsv)
    if isinstance(junc_coord, str):
        if not "-" in junc_coord:
            raise RuntimeError("junc_cood should be '###-###', not %s" % junc_coord)
    elif isinstance(junc_coord, list):
        if isinstance(junc_coord[0], int) and len(junc_coord) == 2:
            junc_coord = "-".join(junc_coord)
        else:
            raise RuntimeError("Bad argument here ...")
    else:
        raise RuntimeError("Unexpected argument")
    lsv_juncs = io_voila_caleb.get_juncs(lsv)
    if not junc_coord in lsv_juncs:
        raise RuntimeError("Junc not in LSV ...")
    indices = index_all(lsv_juncs, junc_coord)
    if len(indices) != 1:
        pdb.set_trace()
    else:
        return indices[0]


def group_lsvs_by_id(lsv_ids):
    """
    Given list of LSV IDs, return dict of Ensembl ids pointing
        at LSVs
    """
    lsv_ids = copy.copy(lsv_ids)
    lsv_ids.sort()
    gene_to_lsvs = dict()
    ensembl_ids = [x.split(":")[0] for x in lsv_ids]
    while len(ensembl_ids) > 0:
        indices = index_all(ensembl_ids, ensembl_ids[0])
        gene_to_lsvs[ensembl_ids[0]] = [lsv_ids[ind] for ind in indices]
        ensembl_ids = ensembl_ids[max(indices) + 1:]
        lsv_ids = lsv_ids[max(indices) + 1:]
    return gene_to_lsvs


def index_all(array, element, nottt=False):
    """Return all indices of array that point to an element.
        Arguments:
            array   = type: list.
            element = object that you wish to get indices for from array.
            Not     = If True, return indices of elements that do *not* match
    """
    if not isinstance(array, list):
        raise Exception("Please only give lists to this function.")
    if nottt:
        matched_indices = [i for i, x in enumerate(array) if x != element]
    else:
        matched_indices = [i for i, x in enumerate(array) if x == element]
    return matched_indices


##################################### GENERAL HELPER #####################################


def percent_through_list(List,
                         By=0.1):
    """
    Given a List (or length of the list), determine which
        indices in the list corrsepond to By percent along the list.

        for example, List = [A,B,C,D,E,F,G,H,I,J] (length = 10)
        if By = 0.2, the dictionary would be:
        {0: 0.0, 2: 20.0, 4: 40.0, 6: 60.0, 8: 80.0, 10: 100.0}

    Use this function to help you determine how many intervals of a given
        percentage along a list loop you are...
    """
    if isinstance(List, list):
        list_n = len(List)
    elif isinstance(List, int):
        list_n = List
        # print "Assuming you gave len(List) instead of the list itself..."
    else:
        raise ValueError("List needs to be a list (or the length of one ...)")
    if not isinstance(By, float) or By <= 0 or By >= 1:
        raise ValueError("By needs to be float >0 and <1")
    # get items at each By percent interval
    by_percent_index = dict()
    for by in io_voila_caleb.calebs_xrange(0, 1, By):
        index = int(math.floor(by * list_n))
        by_percent_index[index] = by * 100.0
    return by_percent_index


def dict_print(obj):
    if type(obj) == dict:
        for k, v in obj.items():
            if hasattr(v, '__iter__'):
                print(k)
                dict_print(v)
            else:
                print('%s : %s' % (k, v))
    elif type(obj) == list:
        for v in obj:
            if hasattr(v, '__iter__'):
                dict_print(v)
            else:
                print(v)
    else:
        print(obj)


def get_inconclusive_cassettes(data, dpsi=0.05, psi=0.05):
    num_d_psis_data = get_num_d_psi(data)
    num_psis_data = get_num_psi(data)
    results = find_binary_lsv_ids(
        num_d_psi=num_d_psis_data,
        this_num_psi=num_psis_data,
        cutoff_d_psi=dpsi,
        cutoff_psi=psi,
        return_other_ids=True,
        return_bi_iis=True)
    binary_ids = results["binary"]
    binary_indices = results["binary_indices"]
    complex_over_ids = results["complex"]
    complex_single_ids = results["single"]
    zero_over_ids = results["non_sig"]
    binary_dict = dict()
    complex_dict = dict()
    single_sig_dict = dict()
    nonsig_dict = dict()

    for comparison in list(data.keys()):
        binary_dict[comparison] = io_voila_caleb.lsv_dict_subset(dictionary=data[comparison],
                                                                 keys=binary_ids,
                                                                 save_LSV_data=True,
                                                                 new_sub_key="binary_indices",
                                                                 new_values=binary_indices)
        complex_dict[comparison] = io_voila_caleb.lsv_dict_subset(dictionary=data[comparison],
                                                                  keys=complex_over_ids,
                                                                  save_LSV_data=True)
        single_sig_dict[comparison] = io_voila_caleb.lsv_dict_subset(dictionary=data[comparison],
                                                                     keys=complex_single_ids,
                                                                     save_LSV_data=True)
        nonsig_dict[comparison] = io_voila_caleb.lsv_dict_subset(dictionary=data[comparison],
                                                                 keys=zero_over_ids,
                                                                 save_LSV_data=True)
    binary_coords = get_binary_coords(binary_dict)
    excl_dict = gen_excl_dict(binary_coords)
    bi_single, bi_cassette, matched_excl_noncassette, bi_complex = get_cassettes(excl_dict)
    all_cassette_ids = list()
    for ii in range(0, len(bi_cassette)):
        cassette = bi_cassette[ii]
        all_cassette_ids.extend(cassette["lsv_ids"])
    return set(all_cassette_ids)


def is_LSV_altSS_weird(LSV):
    """
    Assuming this LSV has a single dPSI change and you want
    to check if all the PSI-util junctions are just alt SS bet
    the ref exon and one other exon.
    """
    psis = np.array(io_voila_caleb.get_psis(LSV))
    chrm = io_voila_caleb.get_chr(LSV)
    strd = io_voila_caleb.get_strand(LSV)
    util_bool = np.max(psis, axis=0) > 0.025
    util_juncs = np.array(io_voila_caleb.get_juncs(LSV))[util_bool].tolist()
    util_juncs = string_to_int_coords(util_juncs)
    ref_ex_int = get_reference_exon(LSV)
    n_ref_ex_int = get_non_reference_exons(LSV)
    n_ref_ex_util_juncs_go_into = list()
    for util_j in util_juncs:
        for side in util_j:
            n_ref_ex_util_juncs_go_into.extend(io_voila_caleb.get_exons_containing_coord(n_ref_ex_int, side))
    unique_exons = set(int_to_string_coords(n_ref_ex_util_juncs_go_into))
    if len(unique_exons) == 1:
        return True
    else:
        return False


def are_used_juncs_sharing_non_ref_coord(Ref_LSV, Junctions):
    """
    Given N junctions (must be 'utilized' by LSV),
        how many of them share the same non-reference-exon coord?

        If N>1 and N of them share non-ref-coord, this is a altSS LSV.

        Else less than N share non-ref-coord, so at least one junction
            is going to a different genomic coordinate. Now, this could
            still be alt SS (esp if they go into same exon), but that's
            outside the scope of this function.
    """
    io_voila_caleb.check_is_lsv(Ref_LSV)
    ref_exon = get_reference_exon(Ref_LSV)
    juncs = [x.split("_")[2] for x in Junctions]
    juncs = string_to_int_coords(juncs)
    non_refs = set()
    poss_differences = len(Junctions)
    if poss_differences == 0:
        pdb.set_trace()
    for junction in juncs:
        non_ref_j = get_non_ref_junc_coord(ref_exon, junction)
        if isinstance(non_ref_j, str):
            continue  # intronic
        non_refs.add(non_ref_j)
    if poss_differences > 1:
        if len(non_refs) == 1:
            return True
        elif len(non_refs) == 0:
            print("likely novel intronic event...")
            pdb.set_trace()
        else:
            # 2+ non-ref coords are different
            return False
    else:  # N= 1, so there is 1 non-ref coord.. duh
        # let us check if other junctions are psi-utilized
        # and show evidence of altSS behavior
        if is_LSV_altSS_weird(Ref_LSV):
            return True  # this means it is a weird LSV with one dPSI
            # and the other junctions are altSS behaving
        return False


def are_used_juncs_sharing_ref_coord(Ref_LSV, Junctions):
    """
    Given N junctions (must be 'utilized' by LSV),
        how many of them share the same reference-exon coord?

        If N>1 and N of them share ref-coord, this is a altSS LSV.

        Else less than N share ref-coord, so at least one junction
            is coming from a different genomic coordinate. Now, this could
            still be alt SS (esp if they go into same exon), but that's
            outside the scope of this function.
    """
    io_voila_caleb.check_is_lsv(Ref_LSV)
    ref_exon = get_reference_exon(Ref_LSV)
    juncs = [x.split("_")[2] for x in Junctions]
    juncs = string_to_int_coords(juncs)
    refs = set()
    poss_differences = len(Junctions)
    if poss_differences == 0:
        pdb.set_trace()
    for junction in juncs:
        ref_j = get_ref_junc_coord(ref_exon, junction)
        if isinstance(ref_j, str):
            continue  # intronic
        refs.add(ref_j)
    if poss_differences > 1:
        if len(refs) == 1:
            return True
        elif len(refs) == 0:
            print("likely novel intronic event...")
            pdb.set_trace()
        else:
            # 2+ ref coords are different
            return False
    else:  # N= 1, so there is 1 ref coord .. duh
        # however:
        # let us check if other junctions are psi-utilized
        # and show evidence of altSS behavior
        if is_LSV_altSS_weird(Ref_LSV):
            return True  # this means it is a weird LSV with one dPSI
        return False


def check_if_lsv_utilizes_intron(data, lsv_id, psi_cutoff=0.05):
    """
    Given ONE LSV, return True if IR PSI is > PSI_cutoff
    """
    io_voila_caleb.check_is_lsv_id(lsv_id)
    LSV = io_voila_caleb.get_lsv(data, lsv_id)
    if LSV["LSV Type"][-1:] != "i":
        return False
    ir_psi = np.array(io_voila_caleb.get_psis(LSV))[:, -1]
    if np.max(ir_psi) > psi_cutoff:
        return True
    return False  # else intron isn't really utilized


def discretize_d_psi(numpy_data,
                     d_psi_thresh=0.1):
    """
    Given a numpy array with rows as junctions, columns as dpsi comparisons,
        modify array so that:
            dpsi >= abs(d_psi_thresh) = 1
            dpsi <= -abs(d_psi_thresh) = -1
            dpsi < abs(d_psi_thresh) = 0
            dpsi > -abs(d_psi_thresh) = 0
    :param numpy_data: numpy array
    :param d_psi_thresh: float threshold. Default 0.1
    :return: nothing
    """
    numpy_data[numpy_data >= abs(d_psi_thresh)] = 1
    numpy_data[numpy_data <= -abs(d_psi_thresh)] = -1
    numpy_data[(numpy_data < abs(d_psi_thresh)) & (numpy_data > -abs(d_psi_thresh))] = 0
    pass


def discretize_prob(numpy_data,
                    prob_thresh=0.0):
    """
    Given a numpy array with rows as junctions, columns as dpsi comparisons,
        modify array so that:
            P(dpsi) >= abs(prob_thresh) = 1
            P(dpsi) < abs(prob_thresh) = 0
    :param numpy_data: numpy array
    :param prob_thresh: float threshold. Default 0.0
    :return: nothing
    """
    numpy_data[numpy_data >= abs(prob_thresh)] = 1
    numpy_data[numpy_data < abs(prob_thresh)] = 0
    pass
