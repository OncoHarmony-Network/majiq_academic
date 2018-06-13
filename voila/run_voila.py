import argparse
import errno
import os
import sys
import time

import voila.constants as constants
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_hdf5 import lsv_id_to_gene_id
from voila.exceptions import VoilaException, CanNotFindVoilaFile
from voila.processes import VoilaPool, VoilaQueue
from voila.utils.utils_voila import create_if_not_exists
from voila.utils.voila_log import voila_log
from voila.view.deltapsi import DeltaPsi
from voila.view.heterogen import Heterogen
from voila.view.psi import Psi
from voila.view.splice_graph import RenderSpliceGraphs


def secs2hms(secs):
    """
    Convert secs into a human readable format.
    :param secs: seconds
    :return: formated time
    """
    m, s = divmod(secs, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def file_versions():
    log = voila_log()
    with SpliceGraph(args.splice_graph) as sg:
        if sg.file_version != constants.SPLICE_GRAPH_FILE_VERSION:
            log.warning('Splice Graph may be out of date.  Verify that you\'re running MAJIQ/Voila with the most '
                        'current versions.')

    if hasattr(args, 'voila_file'):
        for f in args.voila_file:
            with Matrix(f) as m:
                if m.file_version != constants.VOILA_FILE_VERSION:
                    log.warning('Voila file may be out of date.  Verify that you\'re running MAJIQ/Voila with the most '
                                'current versions.')


def gene_names():
    if args.gene_names:
        log = voila_log()
        with SpliceGraph(args.splice_graph) as sg:
            for gene in sg.genes:
                if gene.name in set(args.gene_names):
                    args.gene_ids.append(gene.id)
                    args.gene_names.remove(gene.name)

            if args.gene_names:
                log.warning('Some gene names could not be found in Splice Graph: {}'.format(', '.join(args.gene_names)))

            if not args.gene_ids:
                raise VoilaException('None of the gene names could be converted to gene IDs.')


def gene_ids():
    if args.gene_ids:
        log = voila_log()
        args.gene_ids = list(set(args.gene_ids))
        found_genes = set()
        for f in args.voila_file:
            with Matrix(f) as m:
                for gene_id in args.gene_ids:
                    if any(m.lsv_ids([gene_id])):
                        found_genes.add(gene_id)

        for gene_id in set(args.gene_ids) - found_genes:
            log.warning('Gene ID "{0}" could not be found in a Voila file'.format(gene_id))

        if not found_genes:
            raise VoilaException('None of the gene IDs could be found in a Voila file.')


def lsv_ids():
    if hasattr(args, 'lsv_ids') and args.lsv_ids:
        args.lsv_ids = list(set(args.lsv_ids))
        log = voila_log()
        found_lsvs = set()
        for f in args.voila_file:
            with Matrix(f) as m:
                for lsv_id in args.lsv_ids:
                    if lsv_id in set(m.lsv_ids([lsv_id_to_gene_id(lsv_id)])):
                        found_lsvs.add(lsv_id)

        for lsv_id in set(args.lsv_ids) - found_lsvs:
            log.warning('LSV ID "{0}" could not be found in a Voila file'.format(lsv_id))

        if not found_lsvs:
            raise VoilaException('None of the LSV IDs could be found in a Voila file.')

        args.lsv_ids = list(found_lsvs)


def new_subparser():
    return argparse.ArgumentParser(add_help=False)


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
            raise CanNotFindVoilaFile(value)
        else:
            raise


def check_dir(value):
    """
    check if directory exists.  If not, create it.
    :param value: directory path
    :return: value
    """
    value = os.path.expanduser(value)
    create_if_not_exists(value)
    return value


def check_voila_dir(value):
    return [os.path.join(value, f) for f in os.listdir(value) if f.endswith('.voila')]


def check_file(value):
    """
    Check if file exists.
    :param value: file path
    :return:
    """
    value = os.path.expanduser(value)
    if not os.path.isfile(value):
        raise CanNotFindVoilaFile(value)
    return value


def check_procs(value):
    return min(os.cpu_count(), max(int(value), 1))


parser = argparse.ArgumentParser(description='VOILA is a visualization package '
                                             'for Alternative Local Splicing Events.')
parser.add_argument('-v', action='version', version=constants.VERSION)

# splice graph parser
splice_graph = new_subparser()

splice_graph.add_argument('-s', '--splice-graph', type=check_file, required=True,
                          help='Path to splice graph file.')
splice_graph.add_argument('-o', '--output', type=check_dir, required=True,
                          help='Path for output directory.')
splice_graph.add_argument('--logger', help='Path for log files.')
splice_graph.add_argument('--silent', action='store_true', help='Do not write logs to standard out.')
splice_graph.add_argument('--debug', action='store_true', help='Prints extra output for debugging.')
splice_graph.add_argument('-j', '--nproc', type=check_procs, default=max(int(os.cpu_count() / 2), 1),
                          help='Number of processes used to produce output. Default is half of system processes. ')
splice_graph.add_argument('--gene-names-file', dest='gene_names', type=check_list_file, default=[],
                          help='Location of file that contains a list of common gene names which should remain in '
                               'the results. One name per line.')
splice_graph.add_argument('--gene-names', nargs='*', default=[],
                          help='Common gene names, separated by spaces, which should remain in the results. e.g. '
                               'GENE1 GENE2 ...')
splice_graph.add_argument('--gene-ids-file', dest='gene_ids', type=check_list_file, default=[],
                          help='Location of file that contains a list of gene IDs which should remain in the '
                               'results. One name per line.')
splice_graph.add_argument('--gene-ids', nargs='*', default=[],
                          help='Gene IDs, separated by spaces, which should remain in the results. e.g. GENE_ID1 '
                               'GENE_ID2 ...')

# psi parser
psi_parser = new_subparser()
psi_parser.add_argument('--voila-file', nargs='+', type=check_file,
                        help='Location of majiq\'s voila file.  File should end with ".voila".')
psi_parser.add_argument('--gtf', action='store_true', help='Generate GTF (GFF2) files for LSVs.')
psi_parser.add_argument('--gff', action='store_true', help='Generate GFF3 files for LSVs.')
psi_parser.add_argument('--disable-html', action='store_true', help='Do not write html files.')
psi_parser.add_argument('--disable-tsv', action='store_true',
                        help='Do not generate tab-separated values output file.')
psi_parser.add_argument('--lsv-types-file', type=check_list_file, dest='lsv_types',
                        help='Location of file that contains a list of LSV types which should remain in the results. One '
                             'type per line')
psi_parser.add_argument('--lsv-types', nargs='*', default=[],
                        help='LSV types which should remain in the results')
psi_parser.add_argument('--lsv-ids-file', type=check_list_file, dest='lsv_ids',
                        help='Location of file that contains a list of LSV IDs which should remain in the results. One ID '
                             'per line.')
psi_parser.add_argument('--lsv-ids', nargs='*', default=[],
                        help='LSV IDs, separated by spaces, which should remain in the results. e.g LSV_ID1 LSV_ID2 ...')
psi_parser.add_argument('--disable-db', action='store_true', help='Disables the creation of the database files.')

# deltapsi parser
dpsi_parser = new_subparser()
dpsi_parser.add_argument('--threshold', type=float, default=0.2,
                         help='Filter out LSVs with no junctions predicted to change over a certain value. Even when '
                              'show-all is used this value is still used to calculate the probability in the TSV. The '
                              'default is "0.2".')
dpsi_parser.add_argument('--non-changing-threshold', type=float, default=0.05,
                         help='The default is "0.05".')
dpsi_parser.add_argument('--probability-threshold', type=float, default=None,
                         help='This is off by default.')
dpsi_parser.add_argument('--show-all', action='store_true',
                         help='Show all LSVs including those with no junction with significant change predicted.')

# heterogen parser
het_parser = new_subparser()
het_parser.add_argument('--voila-dir', help='Directory of majiq\'s voila files.', dest='voila_file',
                        type=check_voila_dir)

# subparsers
subparsers = parser.add_subparsers(help='')
subparsers.add_parser('splice-graph', parents=[splice_graph]).set_defaults(func=RenderSpliceGraphs)
subparsers.add_parser('psi', parents=[splice_graph, psi_parser]).set_defaults(func=Psi)
subparsers.add_parser('deltapsi', parents=[splice_graph, psi_parser, dpsi_parser]).set_defaults(func=DeltaPsi)
subparsers.add_parser('heterogen',
                      parents=[splice_graph, psi_parser, dpsi_parser, het_parser]).set_defaults(
    func=Heterogen)

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

args = parser.parse_args()


def main():
    """
    Main function.
    :return: None
    """

    # Time execution time
    start_time = time.time()

    # set up logging
    log_filename = 'voila.log'
    if args.logger:
        log_filename = os.path.join(args.logger, log_filename)
    elif args.output:
        log_filename = os.path.join(args.output, log_filename)
    else:
        log_filename = None

    log = voila_log(filename=log_filename, silent=args.silent, debug=args.debug)
    log.info('Command: {0}'.format(' '.join(sys.argv)))

    log.info('Voila v{}'.format(constants.VERSION))

    VoilaPool(args.nproc)
    VoilaQueue(nprocs=args.nproc)

    try:
        file_versions()
        gene_names()
        gene_ids()
        lsv_ids()

        args.func(args)

        log.info("Voila! Created in: {0}".format(args.output))

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

    VoilaPool().close()


if __name__ == '__main__':
    main()
