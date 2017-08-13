import pdb

from voila.tools import Tool
from voila.utils.voila_log import voila_log
from voila.tools.utils import io_caleb
from voila.tools import remove_dpsi_priors
import pandas as pd
import numpy as np
from voila.tools.utils.percent_through_list import percent_through_list

# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'

LOG = voila_log()


class ThisisLookup(Tool):
    help = 'Given voila tab files and prior removed pickle files, generate matrix files where each row is a' \
           ' lsv\'s junction and each column is a comparison. Matrix files generated will be: E(dpsi), P(E(dpsi)) for' \
           ' each --threshold voila was run at, and prior-removed E(dpsi). If a column is missing data, ' \
           'it will be blank.'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory or file list where voila texts are listed.')
        help_mes = "Output filepath and prefix to save matrices to.)"
        parser.add_argument('outpath',
                            type=str,
                            help=help_mes)
        help_mes = 'Optional: generate a single csv file for all the data? This makes a halfway decent Excel-ready file'
        parser.add_argument('--as_one',
                            action='store_true',
                            default=False,
                            help=help_mes)
        help_mes = 'Optional pattern matching to identify the voila text files'
        parser.add_argument('-p',
                            '--pattern',
                            default="*tsv",
                            type=str,
                            help=help_mes)
        help_mes = "Flag: DO consider IR LSVs (default is to skip)"
        parser.add_argument('--also_ir',
                            action='store_true',
                            help=help_mes,
                            default=False)
        help_mes = "dPSI threshold by which to call junctions as changing"
        parser.add_argument('--dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.2)
        help_mes = "Prob(dPSI) threshold by which to call junctions as changing"
        parser.add_argument('--prob_dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.95)
        help_mes = 'Which comparisons or samples to lookup ID in? Single space or comma separated please.'
        parser.add_argument('--names',
                            '--comparisons',
                            type=str,
                            help=help_mes)
        return parser

    def run(self, args):
        # parse the comparisons argument
        if args.names:
            if "," in args.names or " " in args.names:
                args.names.replace(" ", ",")
                to_lookup = args.names.split(",")
            else:
                to_lookup = [args.names]
        else:
            to_lookup = None
        imported = io_caleb.quick_import(input=args.directory,
                                         pattern=args.pattern,
                                         keep_ir=args.also_ir,
                                         comparisons=to_lookup)
        io_caleb.check_is_ignant(imported, args.dpsi_thresh)
        sig_ids = io_caleb.get_sig_lsv_ids(data=imported,
                                           cutoff_d_psi=args.dpsi_thresh,
                                           prob_d_psi=args.prob_dpsi_thresh,
                                           collapse=True)
        blanked_dict = io_caleb.impute_missing_lsvs(imported,
                                                    impute_with=9)  # imputing with 9s b/c I'll replace this with NA
        io_caleb.quick_import_subset(imported,
                                     sig_ids,
                                     in_place=True)
        # Remove the blanked info for lsvids that were not significantly changing
        for comp in blanked_dict:
            blanked_dict[comp] = set(sig_ids) & blanked_dict[comp]
        flattened = flatten(imported,
                            blank_data=blanked_dict,
                            lsv_ids=sig_ids,
                            return_pandas=not args.as_one)
        if args.as_one:
            dpsi_mat, prb_mat, priormat, row_names = flattened

            all_comparisons = io_caleb.get_comparisons(imported)
            col_order = ["Gene Name", "LSV ID", "Junction"]
            final_frame = pd.DataFrame({keyn: dpsi_mat[keyn] for keyn in col_order})
            final_frame = final_frame[col_order]
            example_lsv = imported[comp][list(sig_ids)[0]]
            prob_header = io_caleb.get_name_of_prob_key(example_lsv)
            for comp in all_comparisons:
                # first 3 cols
                subeaaders = ["dPSI", prob_header, "Confidently non-changing?"]
                header_to_add = [comp + ", " + subheader for subheader in subeaaders]
                col_order.append(header_to_add)
                # data_to_add.extend([dpsi_mat[comp], prb_mat[comp], priormat[comp]])
                thisframe = pd.DataFrame([dpsi_mat[comp], prb_mat[comp], priormat[comp]]).T
                thisframe.columns = header_to_add
                final_frame = pd.concat([final_frame, thisframe], axis=1)
            expanded_header = final_frame.columns.str.split(', ', expand=True).values
            final_frame.columns = pd.MultiIndex.from_tuples([('', x[0]) if pd.isnull(x[1]) else x for x in expanded_header])
            final_frame = final_frame.replace(9, np.NaN)
            final_frame["LSV_JUNC"] = row_names
            final_frame = final_frame.set_index("LSV_JUNC")
            # First row called first, second row called second...
            final_frame.columns.names = ["first", "second"]
            # Multiple dPSI and P(dPSI) by 100
            dpsi_cols = final_frame.columns.get_level_values("second") == "dPSI"
            final_frame.iloc[:, dpsi_cols] = final_frame.iloc[:, dpsi_cols] * 100
            prb_cols = final_frame.columns.get_level_values("second") == prob_header
            final_frame.iloc[:, prb_cols] = final_frame.iloc[:, prb_cols] * 100
            # Make prior removed column just say Yes or be blank
            prior_cols = final_frame.columns.get_level_values("second") == "Confidently non-changing?"
            is_conf_noch = (abs(final_frame.iloc[:, prior_cols]) < 0.05) & (pd.notnull(final_frame.iloc[:, prior_cols]))
            final_frame.iloc[:, prior_cols] = is_conf_noch
            final_frame.iloc[:, prior_cols] = final_frame.iloc[:, prior_cols].replace(False, np.NaN)
            final_frame.iloc[:, prior_cols] = final_frame.iloc[:, prior_cols].replace(1.0, "Yes")

            final_frame.to_csv(args.outpath)
        else:
            dpsi_mat, prb_mat, priormat = flattened
            outp = args.outpath + "dpsi.csv"
            dpsi_mat.to_csv(outp)
            outp = args.outpath + "prob.csv"
            prb_mat.to_csv(outp)
            outp = args.outpath + "priorem.csv"
            priormat.to_csv(outp)


def flatten(imputed_data,
            blank_data,
            lsv_ids,
            return_pandas=False):
    """
    Return a pandas dataframe where each row is a junction and each column is a comparison.
    :param imputed_data:
    :param blank_data: returned from impute_missing_lsvs
    :param lsv_ids:
    :param return_pandas: If True, make a dataframe for each data, Flase, return dicts
    :return: pnadas dataframe
    """
    all_lsvs = io_caleb.get_all_lsv_ids(imputed_data)
    all_comparisons = io_caleb.get_comparisons(imputed_data,
                                               sort=True)
    num_nonchanging = remove_dpsi_priors.get_num_nonchanging(imputed_data,
                                                             blank_info=blank_data,
                                                             impute_with=9,
                                                             as_bools=False)
    all_dpsi_dat = {comp: list() for comp in all_comparisons}
    all_prob_dat = {comp: list() for comp in all_comparisons}
    all_prior_dat = {comp: list() for comp in all_comparisons}
    all_row_names = list()
    lsv_id_col = {"LSV ID": list()}
    junc_col = {"Junction": list()}
    gene_col = {"Gene Name": list()}
    indeces_at_10_percent = percent_through_list(all_lsvs, 0.1)
    i = 0.0
    LOG.info("Flattening %s LSVs ..." % len(all_lsvs))
    for lsvid in all_lsvs:
        if i > 0.0 and i in indeces_at_10_percent:
            LOG.info(str(indeces_at_10_percent[i]) + "% of LSVs flattened...")
        i += 1.0
        gene_name = io_caleb.genename_from_id(imputed_data[all_comparisons[0]], lsvid)
        these_juncs = io_caleb.get_juncs(imputed_data[all_comparisons[0]][lsvid])
        row_names = [("%s_%s" % (lsvid, jj)) for jj in these_juncs]
        [lsv_id_col["LSV ID"].append(lsvid) for jj in these_juncs]
        [junc_col["Junction"].append(jj) for jj in these_juncs]
        [gene_col["Gene Name"].append(gene_name) for jj in these_juncs]
        all_row_names.extend(row_names)
        flat_dpsi = flatten_dpsi(imputed_data, lsvid)
        flat_prob = flatten_prob(imputed_data, lsvid)
        flat_prio = flatten_priors(num_nonchanging, lsvid, all_comparisons)
        [all_dpsi_dat[comp].extend(flat_dpsi[comp]) for comp in all_comparisons]
        [all_prob_dat[comp].extend(flat_prob[comp]) for comp in all_comparisons]
        [all_prior_dat[comp].extend(flat_prio[comp]) for comp in all_comparisons]
    all_dpsi_dat.update(lsv_id_col)
    all_dpsi_dat.update(junc_col)
    all_dpsi_dat.update(gene_col)
    all_prob_dat.update(lsv_id_col)
    all_prob_dat.update(junc_col)
    all_prob_dat.update(gene_col)
    all_prior_dat.update(lsv_id_col)
    all_prior_dat.update(junc_col)
    all_prior_dat.update(gene_col)
    if return_pandas:
        # Make sure order of columns is good
        col_order = ["Gene Name", "LSV ID", "Junction"]
        col_order.extend(all_comparisons)
        dpsi_mat = pd.DataFrame(data=all_dpsi_dat, index=all_row_names)
        dpsi_mat = dpsi_mat.replace(9, np.NaN)
        dpsi_mat = dpsi_mat[col_order]
        prb_mat = pd.DataFrame(data=all_prob_dat, index=all_row_names)
        prb_mat = prb_mat.replace(9, np.NaN)
        prb_mat = prb_mat[col_order]
        priormat = pd.DataFrame(data=all_prior_dat, index=all_row_names)
        priormat = priormat.replace(9, np.NaN)
        priormat = priormat[col_order]
        return dpsi_mat, prb_mat, priormat
    else:
        return all_dpsi_dat, all_prob_dat, all_prior_dat, all_row_names


def flatten_dpsi(imputed_dat, lsvid):
    """
    :param imputed_dat: quick import that was imputed with 9s
    :param lsvid:
    :return: dict of comparison:data to be used for intializing a pandas dataframe
    """
    comparisons = io_caleb.get_comparisons(imputed_dat, sort=True)
    dpsis = [io_caleb.get_dpsis(imputed_dat[comp][lsvid]) for comp in comparisons]
    row_data = {comp: list() for comp in comparisons}
    for dpsi_list, comp in zip(dpsis, comparisons):
        for dpsi in dpsi_list:
            row_data[comp].append(dpsi)
    return row_data


def flatten_prob(imputed_dat, lsvid):
    """

    :param imputed_dat: quick import that was imputed with 9s
    :param lsvid:
    :return: dict of comparison:data to be used for intializing a pandas dataframe
    """
    comparisons = io_caleb.get_comparisons(imputed_dat, sort=True)
    probs = [io_caleb.get_probs(imputed_dat[comp][lsvid]) for comp in comparisons]
    row_data = {comp: list() for comp in comparisons}
    for dpsi_list, comp in zip(probs, comparisons):
        for dpsi in dpsi_list:
            row_data[comp].append(dpsi)
    return row_data


def flatten_priors(nonchanging_np, lsvid, comparisons):
    """

    :param nonchanging_np:
    :param lsvid:
    :param comparisons:
    :return: dict of comparison:data to be used for intializing a pandas dataframe
    """
    # As of now, this array is as follows:
    # rows = junctions, cols = comparisons
    this_nc = nonchanging_np[lsvid]
    as_list = this_nc.T.tolist()
    return {comp: data for comp, data in zip(comparisons, as_list)}
