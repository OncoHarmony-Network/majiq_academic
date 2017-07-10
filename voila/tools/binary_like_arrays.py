from voila.tools import Tool
from voila.tools.utils import io_caleb
from voila.tools import find_binary_lsvs
from voila.tools.utils import random_merge_numpy

import numpy as np
import pandas as pa


# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


class ThisisBinaryLikeArrays(Tool):
    help = 'Given a list of voila dPSI txt files, return matrix of' \
           ' E(PSI), E(dPSI), or P(E(dPS)) from closer, further, or ' \
           'random binary-like LSVs. Rows are LSVs, columns are' \
           ' conditions/comparisons.\n' \
           'Usage:\n' \
           'voila tools binary_like_arrays <path/to/voila/txt/dir> <data_type>'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory where voila texts are.')
        parser.add_argument('data_type',
                            choices={"dpsi", "psi", "prob"},
                            type=str,
                            help='What data should be in the array?')
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
        help_mes = "Choose which junction to get" \
                   " from reference binary-like LSVs."
        parser.add_argument('--junction',
                            type=str,
                            choices={'closer', 'further', 'random'},
                            help=help_mes,
                            default="closer")
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
        help_mes = "Flag: If you don't care whether the binary-like LSVs are reciprocating"
        parser.add_argument('--non_reciprocating',
                            action='store_false',
                            help=help_mes,
                            default=True)
        help_mes = "Output to save file to (include full path)"
        parser.add_argument('-o',
                            '--outfile',
                            type=str,
                            help=help_mes)
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
        imported = io_caleb.quick_import(dir=args.directory,
                                         cutoff_d_psi=0,
                                         cutoff_prob=0,
                                         pattern=args.pattern,
                                         keep_ir=consider_ir)
        if impute:
            blanked_dict = io_caleb.impute_missing_lsvs(data=imported,
                                                        impute_with=0,
                                                        in_place=True,
                                                        warnings=False)
        results_count = find_binary_lsvs.get_binary_lsvs(data=imported,
                                                         method=args.method,
                                                         cutoff_d_psi=args.dpsi_thresh,
                                                         just_lsv_ids=False,
                                                         must_reciprocate=args.non_reciprocating)
        if impute:
            for comp in list(results_count.keys()):
                binary_ids = set(list(results_count[comp].keys()))
                blanked_ids = blanked_dict[comp]
                binary_and_blank_ids = binary_ids & blanked_ids
                blanked_dict[comp] = binary_and_blank_ids
            io_caleb.change_imputed_values(results_count, blanked_dict, new_val=np.NAN)
        the_final_array, lsv_ids, col_names = get_num_array(results_count,
                                                            which_junc=args.junction,
                                                            datatype=args.data_type)
        pandas_df = pa.DataFrame(the_final_array, index=lsv_ids, columns=col_names)
        if args.outfile:
            pandas_df.to_csv(path_or_buf=str(args.outfile),
                             sep="\t")
        else:
            print(pandas_df)


def get_num_array(data, which_junc, datatype):
    """

    :param data: LSV dictionaries returned by find_binary_lsvs.get_binary_lsvs
    :param which_junc: 'closer' or 'further' (wrt reference LSV) or 'random'
    :param datatype: "dpsi", "psi", or "prob"

    row_names = LSV IDs
    col_names = conditions or comparisons (depends on whether PSI or dPSI/Prob)
    :return: (numpy array, row_names, col_names)
    """
    if datatype == "dpsi":
        closer, further, colnames, lsvids = num_d_psi_arrays(data)
    elif datatype == "psi":
        closer, further, colnames, lsvids = num_psi_arrays(data)
    elif datatype == "prob":
        closer, further, colnames, lsvids = num_prob_arrays(data)
    else:
        raise ValueError("Unexpected datatype: %s" % datatype)
    if which_junc == "random":
        rand_merged_array = random_merge_numpy.random_merge(closer, further, preserve_order=True)
        return rand_merged_array, lsvids, colnames
    elif which_junc == "closer":
        return closer, lsvids, colnames
    elif which_junc == "further":
        return further, lsvids, colnames
    else:
        raise ValueError("Unexpected which_junc: %s" % which_junc)


def num_d_psi_arrays(data,
                     return_comparisons=True,  # aka Column names
                     return_lsv_ids=True):  # aka Row names
    """
    Given dictionary of LSV dictionaries that were returned by
    get_binary_LSVs(), return two numpy arrays [n_LSVs, n_conditions]
    where the first array holds dPSI values from the junctions closer to
    the reference exon, and the second array holds dPSI values from
    the junctions further from the reference exon.

    Results returned in the following order:
    lsv_d_psi_arrays_close, lsv_d_psi_arrays_far, Conditions, LSV_IDs
    """
    find_binary_lsvs.check_is_binary_lsv_data(data)
    # dpsi from junction closer to reference
    dpsis_close, conditions = find_binary_lsvs.get_num_d_psi(data,
                                                             return_comparisons=True,
                                                             use_binary_index_info="closer")
    # dpsi from junction further from reference
    dpsis_far, conditions = find_binary_lsvs.get_num_d_psi(data,
                                                           return_comparisons=True,
                                                           use_binary_index_info="further")
    n_cols = len(conditions)
    n_lsvs = len(list(dpsis_close.keys()))
    lsv_d_psi_arrays_close = np.empty([n_lsvs, n_cols])
    lsvs_close = list(dpsis_close.keys())
    lsvs_close.sort()
    for index in range(0, n_lsvs, 1):
        lsv_d_psi_array = dpsis_close[lsvs_close[index]]
        lsv_d_psi_arrays_close[index] = lsv_d_psi_array

    lsv_d_psi_arrays_far = np.empty([n_lsvs, n_cols])
    lsvs_far = list(dpsis_far.keys())
    lsvs_far.sort()
    for index in range(0, n_lsvs, 1):
        lsv_d_psi_array = dpsis_far[lsvs_far[index]]
        lsv_d_psi_arrays_far[index] = lsv_d_psi_array
    # print "The conditions (columns) are: %s"%conditions
    results = list()
    results.append(lsv_d_psi_arrays_close)
    results.append(lsv_d_psi_arrays_far)
    if return_comparisons:
        results.append(conditions)
    if return_lsv_ids:
        results.append(lsvs_close)
    return tuple(results)


def num_psi_arrays(data,
                   return_comparisons=True,  # aka Column names
                   return_lsv_ids=True):  # aka Row names
    """
    Given dictionary of LSV dictionaries that were returned by
    get_binary_LSVs(), return two numpy arrays [n_LSVs, n_conditions]
    where the first array holds PSI values from the junctions closer to
    the reference exon, and the second array holds PSI values from
    the junctions further from the reference exon.

    Results returned in the following order:
    lsv_psi_arrays_close, lsv_psi_arrays_far, Conditions, LSV_IDs
    """
    io_caleb.check_is_binary_lsv_data(data)
    # psi from junction closer to reference
    psis_close, conditions = find_binary_lsvs.get_num_psi(data,
                                                          return_comparisons=True,
                                                          use_binary_index_info="closer")
    # psi from junction further from reference
    psis_far, conditions = find_binary_lsvs.get_num_psi(data,
                                                        return_comparisons=True,
                                                        use_binary_index_info="further")
    n_cols = len(conditions)
    n_lsvs = len(list(psis_close.keys()))
    lsv_psi_arrays_close = np.empty([n_lsvs, n_cols])
    lsvs_close = list(psis_close.keys())
    lsvs_close.sort()
    for index in range(0, n_lsvs, 1):
        lsv_psi_array = psis_close[lsvs_close[index]]
        lsv_psi_arrays_close[index] = lsv_psi_array

    lsv_psi_arrays_far = np.empty([n_lsvs, n_cols])
    lsvs_far = list(psis_far.keys())
    lsvs_far.sort()
    for index in range(0, n_lsvs, 1):
        lsv_psi_array = psis_far[lsvs_far[index]]
        lsv_psi_arrays_far[index] = lsv_psi_array
    # print "The conditions (columns) are: %s"%conditions
    results = list()
    results.append(lsv_psi_arrays_close)
    results.append(lsv_psi_arrays_far)
    if return_comparisons:
        results.append(conditions)
    if return_lsv_ids:
        results.append(lsvs_close)
    return tuple(results)


def num_prob_arrays(data,
                    return_comparisons=True,  # aka Column names
                    return_lsv_ids=True):  # aka Row names
    """
    Given dictionary of LSV dictionaries that were returned by
    get_binary_LSVs(), return two numpy arrays [n_LSVs, n_conditions]
    where the first array holds prob values from the junctions closer to
    the reference exon, and the second array holds prob values from
    the junctions further from the reference exon.

    Results returned in the following order:
    lsv_prob_arrays_close, lsv_prob_arrays_far, Conditions, LSV_IDs
    """
    io_caleb.check_is_binary_lsv_data(data)
    # prob from junction closer to reference
    probs_close, conditions = find_binary_lsvs.get_num_prob(data,
                                                            return_comparisons=True,
                                                            use_binary_index_info="closer")
    # prob from junction further from reference
    probs_far, conditions = find_binary_lsvs.get_num_prob(data,
                                                          return_comparisons=True,
                                                          use_binary_index_info="further")
    n_cols = len(conditions)
    n_lsvs = len(list(probs_close.keys()))
    lsv_prob_arrays_close = np.empty([n_lsvs, n_cols])
    lsvs_close = list(probs_close.keys())
    lsvs_close.sort()
    for index in range(0, n_lsvs, 1):
        lsv_prob_array = probs_close[lsvs_close[index]]
        lsv_prob_arrays_close[index] = lsv_prob_array

    lsv_prob_arrays_far = np.empty([n_lsvs, n_cols])
    lsvs_far = list(probs_far.keys())
    lsvs_far.sort()
    for index in range(0, n_lsvs, 1):
        lsv_prob_array = probs_far[lsvs_far[index]]
        lsv_prob_arrays_far[index] = lsv_prob_array
    # print "The conditions (columns) are: %s"%conditions
    results = list()
    results.append(lsv_prob_arrays_close)
    results.append(lsv_prob_arrays_far)
    if return_comparisons:
        results.append(conditions)
    if return_lsv_ids:
        results.append(lsvs_close)
    return tuple(results)
