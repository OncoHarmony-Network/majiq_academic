import argparse
import errno
import os
import textwrap
import time
from os import path

import voila.constants as constants
from voila.io_voila import Voila
from voila.splice_graphics import SpliceGraph
from voila.utils.utils_voila import create_if_not_exists
from voila.utils.voila_log import voila_log
from voila.view.conditional_table import conditional_table
from voila.view.deltapsi import Deltapsi
from voila.view.lsv_thumbnails import lsv_thumbnails
from voila.view.psi import Psi
from voila.view.splice_graphs import splice_graphs


class VoilaCantFindFile(argparse.ArgumentTypeError):
    def __init__(self, value):
        super(VoilaCantFindFile, self).__init__('cannot find file "{0}"'.format(value))


def secs2hms(secs):
    """
    Convert secs into a human readable format.
    :param secs: seconds
    :return: formated time
    """
    m, s = divmod(secs, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def check_dir(value):
    """
    check if directory exists.  If not, create it.
    :param value: directory path
    :return: value
    """
    create_if_not_exists(value)
    return value


def check_list_file(value):
    """
    Take file, which is a newline separated list of values, and convert it to a list of strings.  Raise error if file
    doesn't exist.
    :param value: file path
    :return: return list of strings
    """
    try:
        with open(value, 'r') as f:
            return [line.strip() for line in f]
    except IOError as e:
        if e.errno == errno.ENOENT:
            raise VoilaCantFindFile(value)
        else:
            raise


def check_file(value):
    """
    Check if file exists.
    :param value: file path
    :return:
    """
    if not path.isfile(value):
        raise VoilaCantFindFile(value)
    return value


def check_splice_graph_file(value):
    check_file(value)
    with SpliceGraph(value, 'r') as sg:
        sg.check_version()
    return value


def check_voila_file(value):
    check_file(value)
    with Voila(value, 'r') as v:
        v.check_version()
    return value


def required_argument(*args, **kwargs):
    parser = args[0]
    required = parser.add_argument_group('required arguments')
    kwargs['required'] = True
    required.add_argument(*args[1:], **kwargs)


def psi_args():
    """
    Psi specific arguments.
    :return: parser
    """

    parser_single = argparse.ArgumentParser(add_help=False)
    return parser_single


def deltapsi_args():
    """
    Deltapsi specific arguments.
    :return: parser
    """
    parser_delta = argparse.ArgumentParser(add_help=False)

    # Probability threshold used to sum the accumulative probability of inclusion/exclusion.
    parser_delta.add_argument('--threshold',
                              type=float,
                              default=0.2,
                              help='Filter out LSVs with no junction predicted to change over a certain value (in '
                                   'percentage).')

    parser_delta.add_argument('--show-all',
                              dest='show_all',
                              action='store_true',
                              default=False,
                              help='Show all LSVs including those with no junction with significant change predicted')

    return parser_delta


def splice_graphs_args():
    """
    Splice graphs specific arguments.
    :return: parser
    """
    parser_splice_graphs = argparse.ArgumentParser(add_help=False)

    parser_splice_graphs.add_argument('splice_graph',
                                      type=check_splice_graph_file,
                                      help='location of majiq\'s splice graph file')

    parser_splice_graphs.add_argument('--limit',
                                      type=int,
                                      default=20,
                                      help='Limit the number of splice graphs shown, default is 20')

    return parser_splice_graphs


def thumbnails_args():
    """
    LSV thumbnail specific arguments.
    :return: parser
    """

    parser_thumbs = argparse.ArgumentParser(add_help=False)
    parser_thumbs.add_argument('--collapsed',
                               action='store_true',
                               default=False,
                               help='Collapsed LSVs thumbnails in the HTML summary.')
    return parser_thumbs


def conditional_table_args():
    """
    Conditional table specific arguments.
    :return: parser
    """

    parser_conditional_table = argparse.ArgumentParser(add_help=False)

    required_argument(
        parser_conditional_table,
        '--cond-pair',
        nargs=2,
        metavar='M1 M2',
        help='condition pair to compare'
    )

    required_argument(
        parser_conditional_table,
        '--sample-files',
        type=check_file,
        nargs='+',
        help='samples Voila output files')

    required_argument(
        parser_conditional_table,
        '--sample-names',
        dest='sample_names',
        nargs='+',
        help='sample names'
    )

    required_argument(
        parser_conditional_table,
        '--pairwise-dir',
        help='root directory where the pairwise delta psi VOILA summaries were created'
    )

    parser_conditional_table.add_argument('--threshold-change',
                                          type=float,
                                          default=0.2,
                                          help='threshold used to filter non-changing LSVs, default is 0.2')
    return parser_conditional_table


def base_args():
    """
    Arguments common to all analysis types.
    :return:  parser
    """

    base_parser = argparse.ArgumentParser(add_help=False)

    required_argument(
        base_parser,
        '-o', '--output',
        dest='output',
        type=check_dir,
        help='path for output directory'
    )

    base_parser.add_argument(
        '-p', '--processes',
        type=int,
        help='Max number of processes used.  Processes will never exceed the number of available processes.'
    )

    base_parser.add_argument(
        '-l', '--logger',
        default=None,
        help='path for log files'
    )

    base_parser.add_argument(
        '-s', '--silent',
        action='store_true',
        help='do not write logs to standard out'
    )
    return base_parser


def html_args():
    html_parser = argparse.ArgumentParser(add_help=False)
    html_parser.add_argument('--no-html',
                             action='store_true',
                             help='do not write html files')

    html_parser.add_argument('--no-tsv',
                             action='store_true',
                             help='do not generate tab-separated values output file')
    return html_parser


def voila_file():
    """
    Arguments needs for analysis types who use the majiq quantifier file.
    :return: parser
    """

    voila_file_parser = argparse.ArgumentParser(add_help=False)

    voila_file_parser.add_argument('voila_file',
                                   type=check_voila_file,
                                   help='location of majiq\'s voila file')

    required_argument(
        voila_file_parser,
        '--splice-graph',
        type=check_splice_graph_file,
        help='path to splice graph file'
    )

    voila_file_parser.add_argument('--gtf',
                                   action='store_true',
                                   help='generate GTF (GFF2) files for LSVs')

    voila_file_parser.add_argument('--gff',
                                   action='store_true',
                                   help='generate GFF3 files for LSVs')

    return voila_file_parser


def gene_search_args():
    """
    Arguments for analysis types who search for specific genes.
    :return: parser
    """

    gene_search_parser = argparse.ArgumentParser(add_help=False)
    gene_search_parser.add_argument('--gene-names-file',
                                    dest='gene_names',
                                    type=check_list_file,
                                    default=[],
                                    help='location of file that contains a list of gene names which should remain in the '
                                         'results, one name per line')

    gene_search_parser.add_argument('--gene-names',
                                    nargs='*',
                                    default=[],
                                    help='gene names which should remain in the results')
    return gene_search_parser


def lsv_type_search_args():
    """
    Arguments for analysis types who search for specific LSV types.
    :return: parser
    """
    lsv_search_parser = argparse.ArgumentParser(add_help=False)
    lsv_search_parser.add_argument('--lsv-types-file',
                                   type=check_list_file,
                                   dest='lsv_types',
                                   help='location of file that contains a list of LSV types which should remain in the '
                                        'results, one type per line')

    lsv_search_parser.add_argument('--lsv-types',
                                   nargs='*',
                                   default=[],
                                   help='LSV types which should remain in the results')

    return lsv_search_parser


def lsv_id_search_args():
    """
    Arguments for analysis types who search for specific LSV IDs.
    :return: parser
    """
    lsv_search_parser = argparse.ArgumentParser(add_help=False)
    lsv_search_parser.add_argument('--lsv-ids-file',
                                   type=check_list_file,
                                   dest='lsv_ids',
                                   help='location of file that contains a list of LSV IDs which should remain in the '
                                        'results, one ID per line')

    lsv_search_parser.add_argument('--lsv-ids',
                                   nargs='*',
                                   default=[],
                                   help='LSV IDs which should remain in the results')

    return lsv_search_parser


def main():
    """
    Main function.
    :return: None
    """

    # Time execution time
    start_time = time.time()

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
            VOILA is a visualization package for Alternative Local Splicing Events.
            -----------------------------------------------------------------------

            '''))

    parser.add_argument('-v', action='version', version=constants.VERSION)

    # Subparser module to agglutinate all subparsers
    subparsers = parser.add_subparsers(dest='type_analysis')
    subparsers.required = True

    majiq_file = voila_file()
    gene_search = gene_search_args()
    lsv_type_search = lsv_type_search_args()
    lsv_id_search = lsv_id_search_args()
    html = html_args()
    base = base_args()

    # Single LSV by Gene(s) of interest
    parser_single = psi_args()
    subparsers.add_parser(constants.ANALYSIS_PSI,
                          help='Single LSV analysis by gene(s) of interest.',
                          parents=[base, html, gene_search, lsv_type_search, lsv_id_search, majiq_file,
                                   parser_single])

    # Delta LSV
    parser_delta = deltapsi_args()
    subparsers.add_parser(constants.ANALYSIS_DELTAPSI,
                          help='Delta LSV analysis by gene(s) of interest.',
                          parents=[base, html, gene_search, lsv_type_search, lsv_id_search, majiq_file,
                                   parser_delta])

    # In-group out-group analysis option
    parser_conditional_table = conditional_table_args()
    subparsers.add_parser(constants.COND_TABLE,
                          help='Generate a HTML table with a list of LSVs changing between conditions in multiple '
                               'samples.',
                          parents=[base, html, gene_search, lsv_type_search, lsv_id_search, parser_conditional_table])

    # Splice graphs generation option
    parser_splicegraphs = splice_graphs_args()
    subparsers.add_parser(constants.SPLICE_GRAPHS,
                          help='Generate only splice graphs.',
                          parents=[base, gene_search, parser_splicegraphs])

    # Thumbnails generation option
    parser_thumbs = thumbnails_args()
    subparsers.add_parser(constants.LSV_THUMBNAILS,
                          help='Generate LSV thumbnails.',
                          parents=[base, lsv_type_search, parser_thumbs])

    # get args
    args = parser.parse_args()

    # set up logging
    log_filename = 'voila.log'
    if args.logger:
        log_filename = os.path.join(args.logger, log_filename)
    log = voila_log(filename=log_filename, silent=args.silent)
    log.info('Starting voila.')

    # set number of processes
    if args.processes and args.processes < constants.PROCESS_COUNT:
        constants.PROCESS_COUNT = args.processes

    if constants.PROCESS_COUNT == 1:
        log.info('Using 1 process.')
    else:
        log.info('Using {0} processes'.format(constants.PROCESS_COUNT))

    # run function for this analysis type
    type_analysis = {
        constants.ANALYSIS_PSI: Psi,
        constants.ANALYSIS_DELTAPSI: Deltapsi,
        constants.SPLICE_GRAPHS: splice_graphs,
        constants.COND_TABLE: conditional_table,
        constants.LSV_THUMBNAILS: lsv_thumbnails
    }

    type_analysis[args.type_analysis](args)

    log.info("Voila! Summaries created in: {0}.".format(args.output))

    # Add elapsed time
    elapsed_str = secs2hms(time.time() - start_time)
    log.info("Execution time: {0}.".format(elapsed_str))


if __name__ == '__main__':
    main()
