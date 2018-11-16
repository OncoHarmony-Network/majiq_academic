import argparse
import errno
import os
import sys
from pathlib import Path

from voila import constants, config
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_hdf5 import lsv_id_to_gene_id
from voila.exceptions import VoilaException, CanNotFindFile
from voila.tsv import Tsv
from voila.utils.utils_voila import create_if_not_exists
from voila.utils.voila_log import voila_log
from voila.view.views import run_service


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
        if hasattr(args, 'voila_file'):
            for f in args.voila_file:
                with Matrix(f) as m:
                    for gene_id in args.gene_ids:
                        if any(m.lsv_ids([gene_id])):
                            found_genes.add(gene_id)

            for gene_id in set(args.gene_ids) - found_genes:
                log.warning('Gene ID "{0}" could not be found in a Voila file'.format(gene_id))

            if not found_genes:
                raise VoilaException('None of the gene IDs could be found in a Voila file.')
        else:
            with SpliceGraph(args.splice_graph) as sg:
                for gene_id in args.gene_ids:
                    if sg.gene(gene_id):
                        found_genes.add(gene_id)

                for gene_id in set(args.gene_ids) - found_genes:
                    log.warning('Gene ID "{0}" could not be found in the Splice Graph file'.format(gene_id))

                if not found_genes:
                    raise VoilaException('None of the gene IDs could be found in the Splice Graph file.')

        args.gene_ids = list(found_genes)


def lsv_ids():
    if hasattr(args, 'lsv_ids') and args.lsv_ids:
        args.lsv_ids = list(set(args.lsv_ids))
        log = voila_log()
        found_lsvs = set()
        for f in args.voila_file:
            with Matrix(f) as m:
                for lsv_id in set(args.lsv_ids) - found_lsvs:
                    if lsv_id in set(m.lsv_ids([lsv_id_to_gene_id(lsv_id)])):
                        found_lsvs.add(lsv_id)

        for lsv_id in set(args.lsv_ids) - found_lsvs:
            log.warning('LSV ID "{0}" could not be found in a Voila file'.format(lsv_id))

        if not found_lsvs:
            raise VoilaException('None of the LSV IDs could be found in a Voila file.')

        args.lsv_ids = list(found_lsvs)
        args.gene_ids = list(lsv_id_to_gene_id(lsv_id) for lsv_id in args.lsv_ids)


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
            raise CanNotFindFile(value)
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
    value = Path(value)

    value = value.expanduser()
    value = value.absolute()

    if value.exists():
        return value
    else:
        raise CanNotFindFile(value)


parser = argparse.ArgumentParser(description='VOILA is a visualization package '
                                             'for Alternative Local Splicing Events.')
parser.add_argument('-v', action='version', version=constants.VERSION)

# tsv parser
tsv_parser = new_subparser()
tsv_parser.add_argument('files', nargs='+', type=check_file,
                        help='List of files or directories which contains the splice graph and voila files.')
tsv_parser.add_argument('-l', '--logger', help='Set log file and location.  There will be no log file if not set.')
tsv_parser.add_argument('--silent', action='store_true', help='Do not write logs to standard out.')
tsv_parser.add_argument('--debug', action='store_true')
tsv_parser.add_argument('-j', '--nproc', type=int, default=min(os.cpu_count(), max(int(os.cpu_count() / 2), 1)),
                        help='Number of processes used to produce output. Default is half of system processes. ')

tsv_parser.add_argument('--threshold', type=float, default=0.2,
                        help='Filter out LSVs with no junctions predicted to change over a certain value. Even when '
                             'show-all is used this value is still used to calculate the probability in the TSV. The '
                             'default is "0.2".')
tsv_parser.add_argument('--non-changing-threshold', type=float, default=0.05,
                        help='The default is "0.05".')

tsv_parser.add_argument('--probability-threshold', type=float, default=None,
                        help='This is off by default.')

tsv_parser.add_argument('--show-all', action='store_true',
                        help='Show all LSVs including those with no junction with significant change predicted.')
tsv_parser.add_argument('-f', '--file-name', required=True, help="Set the TSV file's name and location.")
tsv_parser.add_argument('--lsv-types-file', type=check_list_file, dest='lsv_types',
                        help='Location of file that contains a list of LSV types which should remain in the results. One '
                             'type per line')
tsv_parser.add_argument('--lsv-types', nargs='*', default=[],
                        help='LSV types which should remain in the results')
tsv_parser.add_argument('--lsv-ids-file', type=check_list_file, dest='lsv_ids',
                        help='Location of file that contains a list of LSV IDs which should remain in the results. One ID '
                             'per line.')
tsv_parser.add_argument('--lsv-ids', nargs='*', default=[],
                        help='LSV IDs, separated by spaces, which should remain in the results. e.g LSV_ID1 LSV_ID2 ...')
tsv_parser.add_argument('--gene-names-file', dest='gene_names', type=check_list_file, default=[],
                        help='Location of file that contains a list of common gene names which should remain in '
                             'the results. One name per line.')
tsv_parser.add_argument('--gene-names', nargs='*', default=[],
                        help='Common gene names, separated by spaces, which should remain in the results. e.g. '
                             'GENE1 GENE2 ...')
tsv_parser.add_argument('--gene-ids-file', dest='gene_ids', type=check_list_file, default=[],
                        help='Location of file that contains a list of gene IDs which should remain in the '
                             'results. One name per line.')
tsv_parser.add_argument('--gene-ids', nargs='*', default=[],
                        help='Gene IDs, separated by spaces, which should remain in the results. e.g. GENE_ID1 '
                             'GENE_ID2 ...')

# view parser
view_parser = new_subparser()
view_parser.add_argument('files', nargs='+', type=check_file,
                         help='List of files or directories which contains the splice graph and voila files.')
view_parser.add_argument('--debug', action='store_true')
view_parser.add_argument('-l', '--logger', default=None, help='Path for log files.')
view_parser.add_argument('--silent', action='store_true', help='Do not write logs to standard out.')
view_parser.add_argument('-j', '--nproc', type=int, default=min(os.cpu_count(), max(int(os.cpu_count() / 2), 1)),
                         help='Number of processes used to produce output. Default is half of system processes. ')
view_parser.add_argument('-p', '--port', type=int, default=0,
                         help='Set port to visualize MAJIQ output. Default is a random port.')
view_parser.add_argument('--force-index', action='store_true',
                         help='Create index even if already exists.')

# subparsers
subparsers = parser.add_subparsers(help='')
subparsers.add_parser('tsv', parents=[tsv_parser],
                      help='Generate tsv output for the supplied files.').set_defaults(func=Tsv)
subparsers.add_parser('view', parents=[view_parser],
                      help='Start service to view the visualization for the supplied files.').set_defaults(
    func=run_service)

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

args = parser.parse_args()

log = voila_log(filename=args.logger, silent=args.silent, debug=args.debug)


def main():
    """
    Main function.
    :return: None
    """

    log.info('Command: {0}'.format(' '.join(sys.argv)))

    log.info('Voila v{}'.format(constants.VERSION))

    try:

        config.write(args)
        args.func()

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

    # VoilaPool().close()


if __name__ == '__main__':
    main()
