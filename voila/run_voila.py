import argparse
import os
import sys
import textwrap
import time

import voila.constants as constants
from voila.api import SpliceGraph
from voila.tools import Tools
from voila.utils.exceptions import VoilaException
from voila.utils.voila_log import voila_log
from voila.utils.voila_pool import VoilaPool
from voila.view.conditional_table import ConditionalTable
from voila.view.deltapsi import Deltapsi
from voila.view.heterogen import Heterogen
from voila.view.lsv_thumbnails import LsvThumbnails
from voila.view.psi import Psi
from voila.view.splice_graphs import RenderSpliceGraphs


def secs2hms(secs):
    """
    Convert secs into a human readable format.
    :param secs: seconds
    :return: formated time
    """
    m, s = divmod(secs, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def check_filter(args):
    if hasattr(args, 'gene_names') and args.gene_names:
        if hasattr(args, 'splice_graph') and args.splice_graph:
            with SpliceGraph(args.splice_graph)as sg:
                for gene in sg.genes:
                    if gene.name in args.gene_names:
                        args.gene_ids.append(gene.id)
                        args.gene_names.remove(gene.name)
                        if not args.gene_names:
                            return
        else:
            raise VoilaException(
                'Splice Graph file was not supplied and, therefore, we\'re unable to convert gene names to gene IDs.')


def voila_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
                VOILA is a visualization package for Alternative Local Splicing Events.
                -----------------------------------------------------------------------

                '''))

    parser.add_argument('-v', action='version', version=constants.VERSION)

    subparsers = parser.add_subparsers(dest='type_analysis')
    subparsers.required = True

    subparsers.add_parser(constants.ANALYSIS_PSI,
                          help='Single LSV analysis by gene(s) of interest.',
                          parents=Psi.arg_parents())

    subparsers.add_parser(constants.ANALYSIS_DELTAPSI,
                          help='Delta LSV analysis by gene(s) of interest.',
                          parents=Deltapsi.arg_parents())
    #
    # subparsers.add_parser(constants.COND_TABLE,
    #                       help='Generate a HTML table with a list of LSVs changing between conditions in multiple '
    #                            'samples.',
    #                       parents=ConditionalTable.arg_parents())
    #
    # subparsers.add_parser(constants.SPLICE_GRAPHS,
    #                       help='Generate only splice graphs.',
    #                       parents=RenderSpliceGraphs.arg_parents())
    #
    # subparsers.add_parser(constants.LSV_THUMBNAILS,
    #                       help='Generate LSV thumbnails.',
    #                       parents=LsvThumbnails.arg_parents())
    #
    # subparsers.add_parser(constants.ANALYSIS_HETEROGEN,
    #                       help='Heterogen analysis by gene(s) of interest.',
    #                       parents=Heterogen.arg_parents())

    # parser_tool = subparsers.add_parser(constants.TOOLS,
    #                                     help='Various tools.',
    #                                     parents=Tools.arg_parents())
    # Tools.add_arguments(parser_tool)

    return parser


def main():
    """
    Main function.
    :return: None
    """

    # Time execution time
    start_time = time.time()

    parser = voila_parser()

    # get args
    args = parser.parse_args()

    # set up logging
    log_filename = 'voila.log'
    if hasattr(args, 'logger') and args.logger:
        log_filename = os.path.join(args.logger, log_filename)
    elif hasattr(args, 'output') and args.output:
        log_filename = os.path.join(args.output, log_filename)
    else:
        log_filename = None

    log = voila_log(filename=log_filename, silent=args.silent, debug=args.debug)
    log.info('Command: {0}'.format(' '.join(sys.argv)))

    log.info('Voila {0} v{1}'.format(args.type_analysis, constants.VERSION))

    VoilaPool(args.processes)

    # run function for this analysis type
    type_analysis = {
        constants.ANALYSIS_PSI: Psi,
        constants.ANALYSIS_DELTAPSI: Deltapsi,
        constants.SPLICE_GRAPHS: RenderSpliceGraphs,
        constants.COND_TABLE: ConditionalTable,
        constants.LSV_THUMBNAILS: LsvThumbnails,
        constants.ANALYSIS_HETEROGEN: Heterogen,
        constants.TOOLS: Tools
    }

    try:

        check_filter(args)

        type_analysis[args.type_analysis](args)

        if hasattr(args, 'output'):
            log.info("Voila! Created in: {0}.".format(args.output))

        # Add elapsed time
        elapsed_str = secs2hms(time.time() - start_time)
        log.info("Execution time: {0}.".format(elapsed_str))

    except KeyboardInterrupt:
        log.warning('Voila exiting')

    except VoilaException as ve:
        if args.debug:
            log.exception(ve)
        else:
            log.error(ve)
        exit(1)

    except Exception as e:
        log.exception(e)
        exit(2)


if __name__ == '__main__':
    main()
