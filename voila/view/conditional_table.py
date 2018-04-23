import csv
import os
from collections import defaultdict
from os import path

import numpy as np

from voila.utils.voila_log import voila_log
from voila.view.html import Html


class ConditionalTable(Html):
    def __init__(self, args):
        super().__init__(args)
        conditional_table(args)

    # @classmethod
    # def arg_parents(cls):
    #     parser = cls.get_parser()
    #
    #     cls.required_argument(
    #         parser,
    #         '--cond-pair',
    #         nargs=2,
    #         metavar='M1 M2',
    #         help='Condition pair to compare.'
    #     )
    #
    #     cls.required_argument(
    #         parser,
    #         '--sample-files',
    #         type=cls.check_file,
    #         nargs='+',
    #         help='Samples Voila output files.')
    #
    #     cls.required_argument(
    #         parser,
    #         '--sample-names',
    #         dest='sample_names',
    #         nargs='+',
    #         help='sample names'
    #     )
    #
    #     cls.required_argument(
    #         parser,
    #         '--pairwise-dir',
    #         help='Root directory where the pairwise delta psi VOILA summaries were created.'
    #     )
    #
    #     parser.add_argument('--threshold-change',
    #                         type=float,
    #                         default=0.2,
    #                         help='Threshold used to filter non-changing LSVs.  Default is 0.2.')
    #     parser.add_argument('--best-comparisons',
    #                         type=int,
    #                         help='Filter out all but the best comparisons.  The number remaining is '
    #                              'the user supplied argument.  "Best comparisons" is defined by '
    #                              'comparisons with most agreeing with the least dissagreeing.')
    #
    #     return (
    #         cls.base_args(), cls.html_args(), cls.gene_search_args(), cls.lsv_type_search_args(),
    #         cls.lsv_id_search_args(), cls.output_args(), parser
    #     )


def conditional_table(args):
    """
    Render conditional table output.
    :param args: command line arguments
    :return: None
    """
    lsvs_dict = load_dpsi_tab(args)

    if args.best_comparisons:
        voila_log().info('Only displaying the best {0} comparisons.'.format(args.best_comparisons))
        sorted_lsvs_dict = sorted(lsvs_dict.items(), key=lambda x: (-x[1]['nagree'], x[1]['ndisagree']))
        lsvs_dict = {key: value for key, value in sorted_lsvs_dict[:args.best_comparisons]}

    if not args.no_html:
        render_html(args, lsvs_dict)

    if not args.no_tsv:
        cond_table_tsv(args, lsvs_dict)


def render_html(args, lsvs_dict):
    """
    Render html output.
    :param args: command line arguments
    :param lsvs_dict: dictionary of lsvs
    :return: None
    """
    output_html = "%s_%s_comp_table_%.2f.html" % (args.cond_pair[0], args.cond_pair[1], args.threshold_change)
    sample_names = collect_sample_names(args.sample_names)
    env = Html.get_env()
    sum_template = get_summary_template(args, env)
    table_marks = table_marks_set(len(lsvs_dict))

    log = voila_log()
    log.info("LSVs added to the table: %d" % len(lsvs_dict.keys()))
    log.info("Creating conditional table html...")

    with open(path.join(args.output, output_html), 'w') as voila_output:
        voila_output.write(
            sum_template.render(
                lsvs=lsvs_dict,
                sample_names=sample_names,
                table_marks=table_marks,
                cond_pair=args.cond_pair,
                thres=args.threshold_change
            )
        )

    copy_static(args)


def cond_table_tsv(args, lsvs_dict):
    """
    Create conditional table tsv output.
    :param args: command line arguments
    :param lsvs_dict: dictionary of lsvs
    :return: None
    """

    log = voila_log()
    output_html = get_output_html(args)
    tsv_file = path.join(args.output, output_html.split('.html')[0] + '.tsv')
    sample_names = args.sample_names

    log.info('Creating conditional table TSV...')

    with open(tsv_file, 'w') as csvfile:
        fieldnames = ['Gene', 'LSV ID', '#Disagreeing', '#Agreeing', '#Changing samples'] + sample_names
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()

        for lsv in lsvs_dict:
            row = {
                'Gene': lsvs_dict[lsv]['gene'],
                'LSV ID': lsv,
                '#Disagreeing': lsvs_dict[lsv]['ndisagree'],
                '#Agreeing': lsvs_dict[lsv]['nagree'],
                '#Changing samples': lsvs_dict[lsv]['nchangs']
            }

            for index, sample_name in enumerate(sample_names):
                sample = round(lsvs_dict[lsv]['expecs'][index], 3)
                if sample > -1:
                    row[sample_name] = sample

            writer.writerow(row)


def collect_sample_names(sample_names):
    """
    Create truncated sample names when they're too long.
    :param sample_names: list of sample names
    :return: list
    """
    max_length = 10
    rtn_list = []

    for sample_name in sample_names:
        if len(sample_name) > max_length:
            trunc_sample_name = sample_name[:max_length - 3] + '...'
        else:
            trunc_sample_name = sample_name

        rtn_list.append({'truncated': trunc_sample_name, 'full': sample_name})

    return rtn_list


def load_dpsi_tab(args):
    """
    Load lsv dictionary.
    :param args: command line arguments
    :return: None
    """
    pairwise_dir = args.pairwise_dir
    outdir = args.output
    tab_files_list = args.sample_files
    filter_genes = args.gene_names
    filter_lsvs = args.lsv_ids
    sample_names = args.sample_names
    thres_change = args.threshold_change

    """Load LSV delta psi information from tab-delimited file."""
    lsvs_dict = defaultdict(lambda: defaultdict(lambda: None))
    if pairwise_dir is None:
        pairwise_dir = os.getcwd()

    root_path = None
    path_prefix = '/'.join(['..'] * len(outdir.strip('./').split('/'))) + '/'
    if len(outdir.strip('./')) == 0:
        path_prefix = './'
    # 2-step process:
    #   1. create a data structure finding the most changing junction,
    #   2. select the expected psi from the most changing junction more frequent in the set
    for idx, tab_file in enumerate(tab_files_list):
        with open(tab_file, 'r') as tabf:
            for line in tabf:
                if line.startswith("#"):
                    continue

                fields = line.split()

                if root_path is None:
                    pr, linkk = os.path.split(fields[-1])
                    linkk = linkk.split('#')[0]
                    while not os.path.exists(pairwise_dir + '/' + linkk) and len(pr) > 0:
                        pr, aux = os.path.split(pr)
                        linkk = aux + '/' + linkk

                    if len(pr) == 0:
                        raise Exception('Couldn\'t determine links to delta psi summaries')
                    root_path = pr

                if filter_genes:
                    if fields[0] not in filter_genes and fields[1] not in filter_genes:
                        continue

                if filter_lsvs:
                    if fields[2].upper() not in filter_lsvs:
                        continue

                expecs = [float(aa) for aa in fields[3].split(";")]

                if lsvs_dict[fields[2]]['expecs'] is None:
                    lsvs_dict[fields[2]]['expecs'] = [[]] * len(sample_names)
                    lsvs_dict[fields[2]]['expecs_marks'] = [None] * len(sample_names)
                    lsvs_dict[fields[2]]['links'] = [None] * len(sample_names)
                    lsvs_dict[fields[2]]['njunc'] = [-1] * len(sample_names)

                idx_max = np.argmax([abs(ee) for ee in expecs])

                lsvs_dict[fields[2]]['expecs'][idx] = expecs
                lsvs_dict[fields[2]]['njunc'][idx] = idx_max
                lsvs_dict[fields[2]]['links'][idx] = path_prefix + pairwise_dir + fields[-1].split(root_path)[1]
                lsvs_dict[fields[2]]['gene'] = fields[0]

    for lsv_idx in lsvs_dict.keys():
        if np.max([abs(bb) for ff in lsvs_dict[lsv_idx]['expecs'] for bb in ff]) < thres_change:
            del lsvs_dict[lsv_idx]  # Remove LSVs not passing the changing threshold
            continue

        idx_most_freq = np.argmax(np.bincount(np.array(lsvs_dict[lsv_idx]['njunc'])[
                                                  (np.array(lsvs_dict[lsv_idx]['njunc']) > -1) & np.array(
                                                      [np.any(np.array([abs(fff) for fff in expec]) > thres_change) for
                                                       expec in lsvs_dict[lsv_idx]['expecs']])]))

        # Mark adjusted most changing junction
        lsvs_dict[lsv_idx]['expecs_marks'] = ~np.array(idx_most_freq == lsvs_dict[lsv_idx]['njunc'])

        for idx_exp, expec in enumerate(lsvs_dict[lsv_idx]['expecs']):
            if len(expec) > 0:
                lsvs_dict[lsv_idx]['expecs'][idx_exp] = expec[idx_most_freq]
            else:
                lsvs_dict[lsv_idx]['expecs'][idx_exp] = -1

        lsvs_dict[lsv_idx]['nchangs'] = np.count_nonzero(
            [abs(ee) > thres_change for ee in lsvs_dict[lsv_idx]['expecs'] if ee > -1])

        lsvs_dict[lsv_idx]['njunc'] = idx_most_freq

        exist_expecs = np.array(lsvs_dict[lsv_idx]['expecs'])[(np.array(lsvs_dict[lsv_idx]['expecs']) > -1) & (
                np.array([abs(xx) for xx in lsvs_dict[lsv_idx]['expecs']]) > thres_change)]

        lsvs_dict[lsv_idx]['ndisagree'] = len(exist_expecs) - max(
            (np.count_nonzero(exist_expecs > 0), np.count_nonzero(exist_expecs <= 0)))

        lsvs_dict[lsv_idx]['nagree'] = len(exist_expecs) - min(
            (np.count_nonzero(exist_expecs > 0), np.count_nonzero(exist_expecs <= 0)))

    return lsvs_dict


def get_summary_template(args, env):
    """
    Get summary template for specific type analysis.
    :param args: command line arguments
    :param env: environment variable
    :return: summary template
    """
    template_file_name = args.type_analysis.replace("-", "_") + "_summary_template.html"
    return env.get_template(template_file_name)
