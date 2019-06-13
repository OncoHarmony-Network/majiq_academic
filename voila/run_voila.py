#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path

from voila import constants, config
from voila.exceptions import VoilaException, CanNotFindFile
from voila.tsv import Tsv
from voila.classify import Classify
from voila.voila_log import voila_log
from voila.view.views import run_service


def check_list_file(value):
    """
    Take file, which is a newline separated list of values, and convert it to a list of strings.  Raise error if file
    doesn't exist.
    :param value: file path
    :return: return list of strings
    """

    value = Path(value).expanduser().resolve()

    if value.exists():
        with open(value, 'r') as f:
            return [line.strip() for line in f]
    else:
        raise CanNotFindFile(value)


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

def check_positive(value):
    ivalue = float(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive float value" % value)
    return ivalue


parser = argparse.ArgumentParser(description='VOILA is a visualization package '
                                             'for Alternative Local Splicing Events.')
parser.add_argument('-v', action='version', version=constants.VERSION)

# log parser
log_parser = argparse.ArgumentParser(add_help=False)
log_parser.add_argument('-l', '--logger', help='Set log file and location.  There will be no log file if not set.')
log_parser.add_argument('--silent', action='store_true', help='Do not write logs to standard out.')

# system parser
sys_parser = argparse.ArgumentParser(add_help=False)
sys_parser.add_argument('-j', '--nproc', type=int, default=min(os.cpu_count(), max(int(os.cpu_count() / 2), 1)),
                        help='Number of processes used to produce output. Default is half of system processes. ')
sys_parser.add_argument('--debug', action='store_true')

# webserver parser
webserver_parser = argparse.ArgumentParser(add_help=False)
webserver_parser.add_argument('--web-server', type=str, default='waitress', choices=('waitress', 'gunicorn', 'flask',),
                        help='Web server backend to use. Default is waitress')
webserver_parser.add_argument('--num-web-workers', type=int, default=min(os.cpu_count(), max(int(os.cpu_count() / 2), 1)),
                        help='Number of processes used to handle web I/O (gunicorn workers). '
                             'Only used if the web server is gunicorn. Default is half of system processor count. ')

# tsv parser
tsv_parser = argparse.ArgumentParser(add_help=False)
required_tsv_parser = tsv_parser.add_argument_group('required named arguments')

tsv_parser.add_argument('files', nargs='+', type=check_file,
                        help='List of files or directories which contains the splice graph and voila files.')

required_tsv_parser.add_argument('-f', '--file-name', required=True, help="Output TSV file's name and location.")

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

tsv_parser.add_argument('--lsv-types-file', type=check_list_file, dest='lsv_types',
                        help='Location of file that contains a list of LSV types which should remain in the results. '
                             'One type per line')
tsv_parser.add_argument('--lsv-types', nargs='*', default=[],
                        help='LSV types which should remain in the results')
tsv_parser.add_argument('--lsv-ids-file', type=check_list_file, dest='lsv_ids',
                        help='Location of file that contains a list of LSV IDs which should remain in the results. One '
                             'ID per line.')
tsv_parser.add_argument('--lsv-ids', nargs='*', default=[],
                        help='LSV IDs, separated by spaces, which should remain in the results. e.g LSV_ID1 '
                             'LSV_ID2 ...')
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
view_parser = argparse.ArgumentParser(add_help=False)
view_parser.add_argument('files', nargs='+', type=check_file,
                         help='List of files or directories which contains the splice graph and voila files.')
view_parser.add_argument('-p', '--port', type=int, default=0,
                         help='Set service port. Default is a random.')
view_parser.add_argument('--host', type=str, default='127.0.0.1',
                         help='Set bind address. ex 0.0.0.0 for all interfaces. Default is a 127.0.0.1 (localhost).')
view_parser.add_argument('--force-index', action='store_true',
                         help='Create index even if already exists.')
view_parser.add_argument('--index-file', type=str, default='',
                         help='Location of index file. If specified, will use a separate HDF5 based file for storing '
                              'index data, rather than using input Voila file')
view_parser.add_argument('--strict-indexing', action='store_true',
                         help='When building an index for a study that uses multiple input voila files (such as '
                              'heterogen), verifies that values for each LSV are the same across all input files. '
                              'This protects against accidentally mixing or using inputs from different runs, but '
                              'slows down the indexing process.')
view_parser.add_argument('--splice-graph-only', action='store_true', help=argparse.SUPPRESS)

# classifier parser
classify_parser = argparse.ArgumentParser(add_help=False)
classify_parser.add_argument('files', nargs='+', type=check_file,
                         help='List of files or directories which contains the splice graph and voila files.')
classify_parser.add_argument('--gene-ids', nargs='*', default=[],
                        help='Gene IDs, separated by spaces, which should remain in the results. e.g. GENE_ID1 '
                             'GENE_ID2 ...')
classify_parser.add_argument('--show-all-modules', action='store_true',
                         help='Do not discard modules that are unquantified my Majiq (no LSVs found)')
classify_parser.add_argument('--keep-constitutive', type=int, nargs='?', const=1,
                         help='Do not discard modules with only one junction, implies "--show-all-modules". Turns on '
                              'output of constitutive.tsv and constitutive column in summary output')
classify_parser.add_argument('--putative-multi-gene-regions', action='store_true',
                         help='Only output a single TSV file describing regions found in inputs with complete breaks '
                              'in the gene (no junctions connecting at all). Implies "--keep-constitutive"')
classify_parser.add_argument('--output-complex', action='store_true',
                         help='Complex module data is dumped to all output TSVs, not only summary')
classify_parser.add_argument('--untrimmed-exons', action='store_true',
                         help='Display original Exon coordinates instead of Trimmed coordinates in output TSVs')
classify_parser.add_argument('--decomplexify-psi-threshold', type=float, default=0.01,
                         help='Filter out junctions where PSI is below a certain value (between 0.0 and 1.0). If multiple '
                              'input files are used, only the highest PSI value is used. If 0 (or 0.0) is specified, '
                              'no filtering fill be done. The default is "0.01". (1%%)')
classify_parser.add_argument('--decomplexify-deltapsi-threshold', type=float, default=0.0,
                         help='Filter out junctions where abs(E(dPSI)) is below a certain value (between 0.0 and 1.0). If multiple '
                              'input files are used, only the biggest difference (dPSI) value is used. If 0 (or 0.0) is specified, '
                              'no filtering fill be done. The default is "0.0". ')
classify_parser.add_argument('--decomplexify-reads-threshold', type=int, default=1,
                         help='Filter out junctions where the number of reads is below a certain value (integer). If multiple '
                              'input files are used, only the biggest number of reads is used. The default is "1". ')
classify_parser.add_argument('--min-experiments', default=0.5, type=check_positive,
                        help='Lower threshold for group filters. min_experiments set the minimum number of experiments '
                             'where the different filters check in order to pass an lsv or junction.\n'
                             '\t + <  1 the value is the fraction of the experiments in the group\n'
                             '\t + >= 1 the value is the actual number of experiments. If the number is set to a '
                             'greater number than the size of the group, we use the size instead.\n'
                             '[Default: %(default)s]]')
classify_parser.add_argument('--threshold', type=float, default=0.2,
                        help='Threshold in delta-PSI quantification column. The default is "0.2".')
classify_parser.add_argument('--non-changing-threshold', type=float, default=0.05,
                        help='Threshold in delta-PSI quantification column. The default is "0.05".')
classify_parser.add_argument('--probability-threshold', type=float, default=None,
                        help='This is off by default.')
required_classify_parser = classify_parser.add_argument_group('required named arguments')
required_classify_parser.add_argument('-d', '--directory', required=True, help="All generated TSV files will be dumped in"
                                                                          " this directory")



# subparsers
subparsers = parser.add_subparsers(help='')
subparsers.add_parser('tsv', parents=[tsv_parser, sys_parser, log_parser],
                      help='Generate tsv output for the supplied files.').set_defaults(func=Tsv)
subparsers.add_parser('view', parents=[view_parser, sys_parser, log_parser, webserver_parser],
                      help='Start service to view the visualization for the supplied files.').set_defaults(
    func=run_service)
subparsers.add_parser('classify', parents=[classify_parser, sys_parser, log_parser],
                      help='Classify Splicing events and generate a breakdown of the classification in '
                           'multiple TSV files.').set_defaults(func=Classify)



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


if __name__ == '__main__':
    main()
