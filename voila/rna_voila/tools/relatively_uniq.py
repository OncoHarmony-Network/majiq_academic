import pdb
import numpy as np
from rna_voila.tools import Tool
from rna_voila.tools import remove_dpsi_priors
from rna_voila.tools.find_binary_lsvs import non_redundant_id_dict, counts_per_row
from rna_voila.tools.non_redundant_sets import non_redundant_set
from rna_voila.tools.utils import io_caleb
from rna_voila.tools.utils.percent_through_list import percent_through_list
from rna_voila.voila_log import voila_log

# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'

LOG = voila_log()


class ThisisRelativelyUniq(Tool):
    help = 'Given a directory with voila text files, or a file with a list' \
           ' of voila text file paths: 1) Determine which LSVs are confidently changing, ' \
           ' and are confidently non-changing per comparison. 2) Identify LSVs that are ' \
           'only confidently changing in a given dPSI comparison or multiple comparison(s) ' \
           'with respect to other comparisons where the LSV is confidently non-changing.'

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
        help_mes = "dPSI threshold by which to call junctions as changing"
        parser.add_argument('--dpsi_thresh',
                            '--dpsi',
                            type=float,
                            help=help_mes,
                            default=0.2)
        help_mes = "Prob(dPSI) threshold by which to call junctions as changing"
        parser.add_argument('--prob_dpsi_thresh',
                            '--prob',
                            type=float,
                            help=help_mes,
                            default=0.95)
        help_mes = "Flag: don't consider intron retention LSVs"
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
        help_mes = "Output to save file to (include full path)"
        parser.add_argument('-o',
                            '--outfile',
                            type=str,
                            help=help_mes)
        return parser

    def run(self, args):
        res = relative_uniq_wrapper(directory=args.directory if args.directory else args.file_list,
                                    dpsi_thresh=args.dpsi_thresh,
                                    prob_thresh=args.prob_dpsi_thresh,
                                    pattern=args.pattern)
        if args.outfile:
            with open(args.outfile, "w") as handle:
                handle.write(res)
        else:
            print(res)


def relative_uniq_wrapper(directory,
                          dpsi_thresh,
                          prob_thresh,
                          pattern):
    """"""
    LOG.info("Identifying relatively unique LSVs with E(dPSI) threshold %s, P(E(dPSI) threshold %s ..." % (dpsi_thresh,
                                                                                                           prob_thresh))
    the_lsv_data = io_caleb.quick_import(directory,
                                         pattern=pattern)
    io_caleb.check_is_ignant(the_lsv_data, dpsi_thresh)
    blank_dict = io_caleb.impute_missing_lsvs(data=the_lsv_data, in_place=True, warnings=False)
    num_non_changing = remove_dpsi_priors.get_num_nonchanging(the_lsv_data,
                                                              blank_info=blank_dict)
    num_dpsi = io_caleb.get_num_d_psi(the_lsv_data)
    num_probs, comparisons = io_caleb.get_num_prob(the_lsv_data, return_comparisons=True)
    changing_lsv_ids = io_caleb.get_sig_lsv_ids(data=the_lsv_data,
                                                cutoff_d_psi=dpsi_thresh,
                                                prob_d_psi=prob_thresh,
                                                collapse=True)
    comparisons = np.array(comparisons)
    results = stringent_uniq(num_non_changing,
                             num_dpsi,
                             num_probs,
                             changing_lsv_ids,
                             comparisons,
                             blank_dict,
                             prob_thresh,
                             dpsi_thresh,
                             the_lsv_data)
    printable = ""
    for comp in results["stats"]:
        printable += "Total number of thresh-passing LSVs in %s: %s\n" % (comp,
                                                                          len(io_caleb.get_sig_lsv_ids(
                                                                              the_lsv_data[comp],
                                                                              cutoff_d_psi=dpsi_thresh,
                                                                              prob_d_psi=prob_thresh)))
        printable += "Number of %s-discriminating LSVs with respect to %s\n" % (comp, results["stats"][comp])
        printable += "%s-discriminating LSVs with respect to %s" % (comp, results["comparison"][comp])
    for comp in results["stats"]:
        LOG.info("%s stats: %s" % (comp, results["stats"][comp]))
    return printable


def stringent_uniq(notchanging,
                   dpsidata,
                   probdata,
                   changingdata,
                   comparisons,
                   blankdata,
                   probthresh,
                   dpsithresh,
                   data):
    np.seterr(all='warn')
    import warnings
    warnings.filterwarnings('error')
    results = dict()
    results["comparison"] = dict()
    results["stats"] = dict()
    for lsvid in changingdata:
        # use this variable to subset the numpy arrays such that only
        #  those comparisons (columns) that quantified the LSV are considered
        quantifiable = io_caleb.comparisons_quantifiable(comparisons, blankdata, lsvid)
        if not True in quantifiable:
            # This should NOT happen since these are all significantly changing LSVs...
            raise RuntimeError("Uh oh, %s is not quantifiable in any comparison???" % lsvid)
        # only consider quantifiable comparisons
        dpsidata[lsvid] = dpsidata[lsvid][:, quantifiable]
        probdata[lsvid] = probdata[lsvid][:, quantifiable]
        notchanging[lsvid] = notchanging[lsvid][:, quantifiable]
        ischanging = over_thresh(dpsidata[lsvid],
                                 probdata[lsvid],
                                 dpsithresh,
                                 probthresh)
        # if Trues are only in 1 comparison (column)
        if (ischanging.sum(axis=0) > 0).sum() == 1:
            # TODO : funciton that handles when Trues are only in 1 column
            # Get the dPSI values corresponding to these significanrly changing junctions
            # dpsichanges = dpsidata[lsvid][ischanging]
            # TODO: add function argument specifiying whether only binary-like and/or reciprocal LSVs are considered
            ###
            # These are the row iis corresponding to the junctions that were called changing
            changing_row_iis = (ischanging.sum(axis=1) > 0)
            # These are the col iis corresponding to the cols/comparisons that were called changing
            # There will be one True; the others will be False
            changing_col_iis = (ischanging.sum(axis=0) > 0)
            changing_comparison = comparisons[quantifiable][changing_col_iis]
            # Get all the cols from the above rows in the notchanging array
            not_changing_data_subset = notchanging[lsvid][changing_row_iis, :]
            # Get all the cols from the above rows in the dpsi array
            dpsi_subset = abs(dpsidata[lsvid])[changing_row_iis, :]
            # These are the col iis corresponding to the comparison(s) that didn't pass the changing thresholds
            nonchanging_col_iis = (ischanging.sum(axis=0) == 0)
            # Now, select the columns that correspond to the above indices
            #  This array will be used to identify which comparisons are confidently not changing with respect to the
            # changing comparison
            not_changing_data_subset = not_changing_data_subset[:, nonchanging_col_iis]
            # *all* the nonchangers were unable to be called confidently nonchanging
            if not_changing_data_subset.sum() == 0:
                update_res_dict(results,
                                changing_comparison,
                                np.array(["None"]),
                                lsvid)
                continue
            # Now select the columns that correspond to these indices in the dpsi array too
            dpsi_subset = dpsi_subset[:, nonchanging_col_iis]
            # number of non-changing rows/juncs per col/comparison
            n_nonchange = not_changing_data_subset.sum(axis=0)
            # for each col/comparison, are the # of n_nonchanging equal to the number of juncs changing?
            # if so, then those cols/comparisons are confidently non-changing with respect to the col/comparison that
            # IS changing! woot!
            non_ch_wrt_change = n_nonchange == not_changing_data_subset.shape[0]
            # If the sum below is 0, that means each column had at least one False.
            # Those False(s) are NOT confidently non-changing, meaning they could be changing. They're either
            # (a) dPSI >=0.05 but below user's threshold (taking into account P(E(dPSI)) )
            # or
            # (b) dPSI < 0.05, but when the prior is removed, the dPSI goes to >=0.05
            # If (a), then the junction could really be changing...
            # if (b), well, ... borderline...
            if non_ch_wrt_change.sum() == 0:
                # culprits = dpsi_subset[np.logical_not(not_changing_data_subset)]
                # Check if junctions are case (a) or (b)
                lower_thresh_check = dpsi_subset < 0.05
                # Now, go ahead and re-run some code (TODO: this is horribly inefficient...)
                # number of non-changing (according to simple dPSI check) rows/juncs per col/comparison
                n_nonchange = lower_thresh_check.sum(axis=0)
                # for each col/comparison, are the # of n_nonchanging equal to the number of juncs changing?
                # if so, then those cols/comps are possibly (not confidently) non-changing with respect to the
                #  col/comparison that IS changing ...
                non_ch_wrt_change = n_nonchange == lower_thresh_check.shape[0]
                # Now check again .. this time, if all columns are False, then at least one junction
                # in each comparison has a dPSI >= 0.05 ... case (a) ...
                if non_ch_wrt_change.sum() == 0:
                    continue
            nonchanging_comparisons = comparisons[quantifiable][np.logical_not(changing_col_iis)][non_ch_wrt_change]
            update_res_dict(results,
                            changing_comparison,
                            nonchanging_comparisons,
                            lsvid)
            ### Get more evidence the Falses in not_changing_data_subset are really non-changing
            ### As of now, the not_changing_data_subset array's Trues are confidently non-changing. False could have
            ### border-line cases. If the dpsi array says TODO do I want to do this?
            ### abs(dpsi_subset) <= 0.05
        else:
            continue
            # else Trues not just in 1 column
            # TODO: something
    return results


def update_res_dict(results,
                    changing_comp,
                    nonchanging_comps,
                    lsv_id):
    """
    Update results dictionary to include stats and lsv ids
    :param results:
    :param changing_comp:
    :param nonchanging_comps:
    """
    changing_comp = list(changing_comp)[0]
    nonchanging_comps = list(nonchanging_comps)
    nonchanging_comps.sort()
    nonchanging_comps = " + ".join(nonchanging_comps)
    if changing_comp not in results["comparison"]:
        results["comparison"][changing_comp] = dict()
        results["stats"][changing_comp] = dict()
    if nonchanging_comps not in results["comparison"][changing_comp]:
        results["comparison"][changing_comp][nonchanging_comps] = list()
        results["stats"][changing_comp][nonchanging_comps] = 0
    results["comparison"][changing_comp][nonchanging_comps].append(lsv_id)
    results["stats"][changing_comp][nonchanging_comps] += 1


def over_thresh(num_dpsi, num_prob, dpsi_thresh, prob_thresh):
    """
    :param num_dpsi: rows = junctions, cols = comparisons
    :param num_prob: rows = junctions, cols = comparisons
    :param dpsi_thresh: abs(num_dpsi) >= abs(dpsi_thresh)
    :param prob_thresh: num_prob >= prob_thresh
    :return: numpy array of True/False
    """
    bool_array = (abs(num_dpsi) >= abs(dpsi_thresh)) & (num_prob >= prob_thresh)
    return bool_array


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
        the results pass the ole' "Yeah that checks out" intuition view.

    Arguments:
        data : Quick Import
        threshold (default_view 0.05) : dPSI w/in +/-[threshold] of 0 are converted to 0
        sig_dpsi_thresh (default_view 0.2) : dPSI thresh for ID'ing sig LSVs
        sig_prob_thresh (default_view 0.95) : Prob(dPSI) thresh for ID'ing sig LSVs
        evaluate_introns : True/False, should intronic LSVs be considered?
        comparisons (default_view False) : return list of comparisons corresponding to columns?
        stringent_prob (default_view 0) : convert dPSI to 0 if P(dPSI) < stringent_prob_thresh
            default_view 0 means no dPSI will be converted to 0 by this argument
            0.95 is very stringent, and will vastly reduce power to detect unique changes
                however, these unique changes will be highly supported by the RNA seq data!
        must_reciprocate (default_view True) : must the uniquely called LSV dPSIs sum to 0?
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

     --(4)--
     --Detect groups of comparisons that agree/disagree (TBD)
    if:
     - at least one row is assigned more than one column name
     - all rows agree on column name assignments
    Then that LSV is assigned to a combination of the comparisons
    eg: Comparison1_Comparison2 <- names will be sorted before joined

    --(5)--
    --Group LSVs together if they belong to the same non-redundant set

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
    LOG.info("Building non_red id dictionary ...")
    id_dict = non_redundant_id_dict(non_red_sets,
                                    cutoff_dpsi=0.1,
                                    cutoff_psi=1)
    LOG.info("Getting all Prob(dPSI) data ...")
    all_prob_data, comp_names = io_caleb.get_num_prob(data, True)
    sig_subset = io_caleb.subset_significant(data,
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
            these_ids = io_caleb.get_lsv_ids(sig_subset[c])
            sig_subset_ids[c] = these_ids
            all_sig.extend(these_ids)
        comparisons_ii = [comp_names.index(x) for x in comparisons]
        all_sig = list(set(all_sig))
    else:
        for c in data.keys():
            sig_subset_ids[c] = io_caleb.get_lsv_ids(sig_subset[c])
        all_sig = io_caleb.get_all_lsv_ids(sig_subset)
        comparisons_ii = [comp_names.index(x) for x in data.keys()]
    thelsvs = list(all_num_dpsis_data.keys())
    indeces_at_10_percent = percent_through_list(thelsvs, 0.1)
    LOG.info("Assigning %s LSVs to groups ... " % (len(thelsvs)))
    i = 0.0
    for lsv_id in thelsvs:
        if i > 0.0 and i in indeces_at_10_percent:
            LOG.info(str(indeces_at_10_percent[i]) + "% of LSVs processed...")
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
                    LOG.info("unsure0")
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
                LOG.info("unsure4")
                pdb.set_trace()

    LOG.info("Finished analyzing LSVs!")
    all_uniques = list(set(all_uniques))
    comp_names = list(singly_unique_sig.keys())
    comp_names.sort()
    nonrednetworks = dict()
    nonrednetworks["all_non_red_sets"] = non_red_sets
    nonrednetworks["unique_nonredsets"] = dict()
    nonrednetworks["cococombos"] = {x: [] for x in cococombo_breaker.keys()}
    nonrednetworks["opposites"] = {x: [] for x in opposite_dict.keys()}
    all_non_red_sets = list(non_red_sets['singles'].values())
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
    LOG.info(summary_text)
    if unblank_the_data:
        io_caleb.unimpute_lsv_data(data, blanked_lsvs_dict)
    results = {"singly_unique": singly_unique_sig,
               "cococombos": cococombo_breaker,
               "opposites": opposite_dict,
               "nonredsets": nonrednetworks,
               "summary": summary_text}
    return results
