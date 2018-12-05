from voila.tools import Tool, non_redundant_sets
from voila.tools.find_binary_lsvs import check_is_binary_lsv_data
from voila.tools.utils import io_caleb
from voila.tools import find_binary_lsvs
from voila.tools.utils import random_merge_numpy
from voila.tools.utils.percent_through_list import percent_through_list
from voila.voila_log import voila_log
import numpy as np
import pandas as pa
import copy
import pickle as pkl

# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'

LOG = voila_log()


class ThisisBinaryLikeArrays(Tool):
    help = 'Given a list of voila dPSI txt files, return matrix of' \
           ' E(PSI), E(dPSI), or P(E(dPS)) from closer, further, or ' \
           'random binary-like LSVs. Rows are LSVs, columns are' \
           ' dPSI comparisons.\n' \
           'Usage:\n' \
           'voila tools binary_like_arrays <path/to/voila/txt/dir> <data_type>'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory where voila texts are OR file with list of text files.')
        parser.add_argument('data_type',
                            choices={"dpsi", "psi", "prob"},
                            type=str,
                            help='What data should be in the array?')
        help_mes = "dPSI threshold by which to call junctions as changing"
        parser.add_argument('--dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0)
        help_mes = "Prob(dPSI) threshold by which to call junctions as changing"
        parser.add_argument('--prob_dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.0)
        help_mes = "PSI threshold by which to call junctions as changing. Default 1 means it isn't used"
        parser.add_argument('--psi_thresh',
                            type=float,
                            help=help_mes,
                            default=1)
        help_mes = "Choose which junction to get" \
                   " from reference binary-like LSVs."
        parser.add_argument('--junction',
                            type=str,
                            choices={'closer', 'further', 'random'},
                            help=help_mes,
                            default="closer")
        help_mes = "Flag: also include IR LSVs (default_view: don't import them)"
        parser.add_argument('--also_ir',
                            action='store_true',
                            help=help_mes,
                            default=False)
        help_mes = "Flag: don't use IR LSVs for identifying non-redundant sets"
        parser.add_argument('--no_non_red_ir',
                            action='store_false',
                            help=help_mes,
                            default=True)
        help_mes = 'Optional pattern matching to identify the voila text files'
        parser.add_argument('-p',
                            '--pattern',
                            default="*tsv",
                            type=str,
                            help=help_mes)
        help_mes = 'Remove LSVs in corresponding genes from results before analysis. ' \
                   'This arg should be a file path to a pickle file. dict[comparisons] -> Gene IDs'
        parser.add_argument('--remove_these_genes',
                            type=str,
                            help=help_mes)
        help_mes = 'Method: \'count\' or \'sum_to_95\''
        parser.add_argument('-m',
                            '--method',
                            choices={"sum_to_95", "count"},
                            default="sum_to_95",
                            type=str,
                            help=help_mes)
        help_mes = "Flag: If you think the binary-like LSVs should reciprocate\n" \
                   "Note: this is broken at the moment, so don't use it.."
        parser.add_argument('--must_reciprocate',
                            action='store_true',
                            help=help_mes,
                            default=False)
        help_mes = "Output to save file to (include full path)"
        parser.add_argument('-o',
                            '--outfile',
                            type=str,
                            help=help_mes)
        return parser

    def run(self, args):
        if not args.no_non_red_ir:
            raise RuntimeError("Why not use IR LSVs for non-red set identification??"
                               " Caleb really doesn't think you should do this, "
                               "so he hasn't implemented it yet.")
        if args.remove_these_genes:
            LOG.info("User will ignore LSVs using data from %s" % args.remove_these_genes)
            rem_dict = pkl.load(open(args.remove_these_genes, 'rb'), encoding='bytes')
        else:
            rem_dict = None
        imported = io_caleb.quick_import(input=args.directory,
                                         cutoff_d_psi=0,
                                         cutoff_prob=0,
                                         pattern=args.pattern,
                                         keep_ir=args.also_ir,
                                         remove_these_genes=rem_dict)
        io_caleb.check_is_ignant(imported, args.dpsi_thresh)
        nrset, blanked_dict = non_redundant_sets.non_redundant_set(data=imported,
                                                                   cutoff_dpsi=args.dpsi_thresh,
                                                                   cutoff_psi=args.psi_thresh,
                                                                   save_blanked_structure=True,
                                                                   bi_method=args.method)
        sig_ids = io_caleb.get_sig_lsv_ids(imported,
                                           cutoff_d_psi=args.dpsi_thresh,
                                           prob_d_psi=args.prob_dpsi_thresh,
                                           collapse=True)
        results_count = find_binary_lsvs.get_binary_lsvs(data=imported,
                                                         method=args.method,
                                                         cutoff_d_psi=None if args.dpsi_thresh == 0 else args.dpsi_thresh,
                                                         cutoff_psi=None if args.psi_thresh == 1 else args.psi_thresh,
                                                         just_lsv_ids=False,
                                                         must_reciprocate=args.must_reciprocate)

        all_binary_ids = set()
        for comp in list(results_count.keys()):
            binary_ids = set(list(results_count[comp].keys()))
            all_binary_ids = all_binary_ids | binary_ids
            blanked_ids = blanked_dict[comp]
            binary_and_blank_ids = binary_ids & blanked_ids
            blanked_dict[comp] = binary_and_blank_ids
        sig_and_binary = list(set(sig_ids) & all_binary_ids)
        non_red_lsv_ids = remove_redundants(imported, nrset, sig_and_binary)
        io_caleb.change_imputed_values(results_count, blanked_dict, new_val=np.NAN)
        the_final_array, bi_lsv_ids, col_names = get_num_array(results_count,
                                                            which_junc=args.junction,
                                                            datatype=args.data_type)
        pandas_df = pa.DataFrame(the_final_array, index=bi_lsv_ids, columns=col_names)
        # Remove redundant rows
        pandas_df = pandas_df.loc[pandas_df.index.isin(non_red_lsv_ids)]
        if args.outfile:
            pandas_df.to_csv(path_or_buf=str(args.outfile),
                             sep="\t")
        else:
            print_full(pandas_df)


def remove_redundants(imputed_data, nrset, lsv_ids):
    """
    Given Imputed data, non-redundant set results, and a list of LSV IDs,
        determine which of your LSV IDs are redundant, and only keep those that
        have the biggest dPSI value.
    :param imputed_data: blanked quick import
    :param nrset: non_red_set
    :param lsv_ids: LSV IDs that have redundancy
    :return: subset of lsv_ids such that none overlap
    """
    non_red_ids = list()
    iis_at_1_percent = percent_through_list(lsv_ids, 0.01)
    i = 0.0
    LOG.info("Removing redundants IDs from list of %s LSV IDs" % len(lsv_ids))
    while len(lsv_ids) > 0:
        if i > 0.0 and i in iis_at_1_percent:
            LOG.info(str(iis_at_1_percent[i]) + "%% of overlapping sets processed (%s LSVs)" % i)
        this_id = lsv_ids.pop()
        # Determine which, if any, LSVs overlap with this_id in your list of lsv_ids
        partners = non_redundant_sets.find_set_partners(nrset, this_id, lsv_ids)
        if partners in ["isolated", "no_sig_partner", "no_sig_partners"]:
            non_red_ids.append(this_id)
            i += 1.0
            continue
        # Ok, this_id has at least 1 partner.
        overlapping_ids = copy.copy(partners)
        overlapping_ids.append(this_id)
        # Only care about testing incoming LSVs
        overlapping_ids = list(set(overlapping_ids) & set(lsv_ids))
        most_changing_id = non_redundant_sets.most_changing_lsv(imputed_data, overlapping_ids)
        non_red_ids.append(most_changing_id)
        for lid in list(set(overlapping_ids) - {this_id}):
            lsv_ids.remove(lid)
        i += 1.0
    return non_red_ids


def print_full(x):
    pa.set_option('display.max_rows', len(x))
    print(x)
    pa.reset_option('display.max_rows')


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
    dpsis_close, conditions = io_caleb.get_num_d_psi(data,
                                                     return_comparisons=True,
                                                     use_binary_index_info="closer")
    # dpsi from junction further from reference
    dpsis_far, conditions = io_caleb.get_num_d_psi(data,
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
    check_is_binary_lsv_data(data)
    # psi from junction closer to reference
    psis_close, conditions = io_caleb.get_num_psi(data,
                                                  return_comparisons=True,
                                                  use_binary_index_info="closer")
    # psi from junction further from reference
    psis_far, conditions = io_caleb.get_num_psi(data,
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
    check_is_binary_lsv_data(data)
    # prob from junction closer to reference
    probs_close, conditions = io_caleb.get_num_prob(data,
                                                    return_comparisons=True,
                                                    use_binary_index_info="closer")
    # prob from junction further from reference
    probs_far, conditions = io_caleb.get_num_prob(data,
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
