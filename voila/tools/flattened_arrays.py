import pdb

from voila.tools import Tool
from voila.tools.utils import io_caleb
from voila.tools import remove_dpsi_priors
import pandas as pd
import numpy as np


# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


class ThisisLookup(Tool):
    help = 'Directory with voila tab files and prior removed pickle files, generate matrix files where each row is a' \
           ' lsv\'s junction and each column is a comparison. Matrix files generated will be: E(dpsi), P(E(dpsi)) for' \
           ' each --threshold voila was run at, and prior-removed E(dpsi). If a column is missing data, it will say NA.'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory or file list where voila texts are listed.')
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
        help_mes = "Flag: don't consider IR LSVs"
        parser.add_argument('--no_ir',
                            action='store_true',
                            help=help_mes,
                            default=False)
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
            dont_remove_dups = False
        else:
            to_lookup = None
            dont_remove_dups=True
        imported = io_caleb.quick_import(input=args.directory,
                                         pattern=args.pattern,
                                         keep_ir=True,
                                         comparisons=to_lookup)
        blanked_dict = io_caleb.impute_missing_lsvs(imported,
                                                    impute_with=9)  # imputing with 9s b/c I'll replace this with NA
        flat_info = flatten(imported, blank_data=blanked_dict)



def flatten(imputed_data,
            blank_data):
    """
    Return a pandas dataframe where each row is a junction and each column is a comparison.
    :param imputed_data:
    :param blank_data: returned from impute_missing_lsvs
    :return: pnadas dataframe
    """
    all_lsvs = io_caleb.get_all_lsv_ids(imputed_data)
    all_comparisons = io_caleb.get_comparisons(imputed_data,
                                               sort=True)
    num_nonchanging = remove_dpsi_priors.get_num_nonchanging(imputed_data,
                                                             blank_info=blank_data,
                                                             impute_with=9,
                                                             as_bools=False)
    all_dpsi_dat = dict()
    all_prob_dat = dict()
    all_prior_dat = dict()
    for lsvid in all_lsvs:
        these_juncs = io_caleb.get_juncs(imputed_data[all_comparisons[0]][lsvid])
        row_names = [("%s_%s" % (lsvid, jj)) for jj in these_juncs]
        flat_dpsi = flatten_dpsi(imputed_data, lsvid)
        flat_prob = flatten_prob(imputed_data, lsvid)
        flat_prio = flatten_priors(num_nonchanging, lsvid, all_comparisons)
        all_dpsi_dat.update(flat_dpsi)
        all_prob_dat.update(flat_prob)
        all_prior_dat.update(flat_prio)

    res = pd.DataFrame(data=row_data, index=row_names)
    res = res.replace(9, np.NaN)

def flatten_dpsi(imputed_dat, lsvid):
    """
    :param imputed_dat: quick import that was imputed with 9s
    :param lsvid:
    :return: dict of comparison:data to be used for intializing a pandas dataframe
    """
    comparisons = io_caleb.get_comparisons(imputed_dat, sort=True)
    dpsis = [io_caleb.get_dpsis(imputed_dat[comp][lsvid]) for comp in comparisons]
    row_data = {comp: list() for comp in comparisons}
    for dpsi_list,comp in zip(dpsis,comparisons):
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
    for dpsi_list,comp in zip(probs,comparisons):
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
    as_list = this_nc.tolist()
    res = {comp}
    pdb.set_trace()
