#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path

from rna_voila import constants, config
from rna_voila.exceptions import VoilaException, CanNotFindFile
from rna_voila.tsv import Tsv
from rna_voila.voila_log import voila_log
from rna_voila.view.views import run_service
from rna_voila.classify import Classify
from rna_voila.filter import Filter
from rna_voila.splitter import splitter, recombine

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
tsv_parser.add_argument('--show-read-counts', action='store_true', help=argparse.SUPPRESS)

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
view_parser.add_argument('--skip-type-indexing', action='store_true',
                         help='Skips creating index for lsv type data (alt3, alt5, binary, etc). These filters will'
                              ' no longer function, but there will be a significant indexing speedup')
view_parser.add_argument('--strict-indexing', action='store_true',
                         help='When building an index for a study that uses multiple input voila files (such as '
                              'heterogen), verifies that values for each LSV are the same across all input files. '
                              'This protects against accidentally mixing or using inputs from different runs, but '
                              'slows down the indexing process.')
view_parser.add_argument('--splice-graph-only', action='store_true', help=argparse.SUPPRESS)
view_parser.add_argument('--enable-passcode', action='store_true',
                         help='Disallow access to the viewer unless the special link is used to start the session'
                              ' (provides some security against port scanners accessing your app')

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
classify_parser.add_argument('--keep-no-lsvs', action='store_true',
                         help='If there are no LSVs attached to a specific junction, retain the junction instead of removing it')
classify_parser.add_argument('--putative-multi-gene-regions', action='store_true',
                         help='Only output a single TSV file describing regions found in inputs with complete breaks '
                              'in the gene (no junctions connecting at all). Implies "--keep-constitutive"')
classify_parser.add_argument('--output-complex', action='store_true',
                         help='Complex module data is dumped to all output TSVs, not only summary')
classify_parser.add_argument('--untrimmed-exons', action='store_true',
                         help='Display original Exon coordinates instead of Trimmed coordinates in output TSVs')
classify_parser.add_argument('--decomplexify-psi-threshold', type=float, default=None,
                         help='Filter out junctions where PSI is below a certain value (between 0.0 and 1.0). If multiple '
                              'input files are used, only the highest PSI value is used. If 0 (or 0.0) is specified, '
                              'no filtering fill be done. The default is "0.01". (1%%)')
classify_parser.add_argument('--decomplexify-deltapsi-threshold', type=float, default=None,
                         help='Filter out junctions where abs(E(dPSI)) is below a certain value (between 0.0 and 1.0). If multiple '
                              'input files are used, only the biggest difference (dPSI) value is used. If 0 (or 0.0) is specified, '
                              'no filtering fill be done. The default is "0.0". ')
classify_parser.add_argument('--decomplexify-reads-threshold', type=int, default=1,
                         help='Filter out junctions where the number of reads is below a certain value (integer). If multiple '
                              'input files are used, only the biggest number of reads is used. The default is "%(default)s". ')

classify_parser.add_argument('--non-changing-pvalue-threshold', type=float, default=0.05,
                        help='For determining non-changing with HET inputs. Minimum p-value for which an LSV/junction'
                             ' can return true. Uses minimum p-value from all tests provided. The default is "%(default)s".')
classify_parser.add_argument('--non-changing-within-group-IQR', type=float, default=0.1,
                        help='For determining non-changing with HET inputs. Maximum IQR within a group for which an '
                             'LSV/junction can return true. The default is "%(default)s".')
classify_parser.add_argument('--non-changing-between-group-dpsi', type=float, default=0.05,
                        help='For determining non-changing with HET inputs. Maximum absolute difference in median '
                             'values of PSI for which an LSV/junction can return true. The default is "%(default)s".')
classify_parser.add_argument('--changing-pvalue-threshold', type=float, default=0.05,
                        help='For determining changing with HET inputs. Maximum p-value for which an LSV/junction'
                             ' can return true. Uses maximum p-value from all tests provided. The default is "%(default)s".')
classify_parser.add_argument('--changing-between-group-dpsi', type=float, default=0.2,
                        help='For determining changing with HET or delta-PSI inputs. For HET, minimum absolute difference in median '
                             'values of PSI for which an LSV/junction can return true. For delta-PSI, min(E(dPSI)).'
                             ' The default is "%(default)s".')
classify_parser.add_argument('--changing-between-group-dpsi-secondary', type=float, default=0.1,
                        help='Set the secondary changing event definition. In order to be considered "changing", and junction in an event must'
                             ' meet the other changing definitions, and ALL junctions in an event must meet this condition (DPSI value'
                             ' of the junction >= this value)'
                             ' The default is "%(default)s".')
classify_parser.add_argument('--non-changing-threshold', type=float, default=0.05,
                        help='Threshold in delta-PSI quantification column. The default is "%(default)s".')
classify_parser.add_argument('--probability-changing-threshold', type=float, default=0.95,
                        help='The default is "%(default)s"')
classify_parser.add_argument('--probability-non-changing-threshold', type=float, default=0.95,
                        help='The default is "%(default)s"')

classify_parser.add_argument('--changing', action='store_true',
                         help='In general, find classifications for events which are changing (between multiple analysis). Requires at least one deltapsi or het voila file')
classify_parser.add_argument('--heatmap-selection', choices=['shortest_junction', 'max_abs_dpsi'],
                         help='For the classifier output "heatmap", the quantification values may be derived from either the shortest junction in the module (default), '
                              'or optionally, if a het or dpsi file is provided, from the junction with the maximum dpsi value')
classify_parser.add_argument('--debug-num-genes', type=int,
                         help='Modulize only n many genes, mainly used for debugging')
classify_parser.add_argument('--overwrite', action='store_true',
                         help='If there are files inside of specified --directory, delete them and run classifier anyway')
classify_parser.add_argument('--enabled-outputs', type=str,
                         help='By default, classifier makes summary output only. However, additional outputs may be enabled. Please check the '
                              'documentation for more information on them. Specify a comma separated list with one or more of the following options: '
                              'summary,events,junctions,heatmap. You may also specify "all" to enable all outputs. If ran with --output-training-data,'
                              'accepted outputs are instead "junctions", "exons" and/or "matrices" (and "all" still works too). In this case default is "all".')
required_classify_parser = classify_parser.add_argument_group('required named arguments')
required_classify_parser.add_argument('-d', '--directory', required=True, help="All generated TSV files will be dumped in"
                                                                          " this directory")

# filter parser
filter_parser = argparse.ArgumentParser(add_help=False)
filter_parser.add_argument('files', nargs='+', type=check_file,
                         help='List of files or directories which contains the splice graph and voila files.')
filter_parser.add_argument('--overwrite', action='store_true',
                         help='If the output filename already exists in the destination directory, overwrite'
                              'it instead of skipping with a warning')
filter_parser.add_argument('--voila-files-only', action='store_true',
                         help='Only filter the voila files, not the splicegraph')
filter_parser.add_argument('--splice-graph-only', action='store_true',
                         help='Only filter the splicegraph, not the voila files')
filter_parser.add_argument('--gene-ids', nargs='*', default=[],
                        help='Gene IDs, separated by spaces, which should remain in the results. e.g. GENE_ID1 '
                             'GENE_ID2 ...')
filter_parser.add_argument('--gene-ids-file', type=check_file,
                        help='Specify Gene IDs to retain in a file instead of on the command line. One Gene ID per line')
filter_parser.add_argument('--lsv-ids', nargs='*', default=[],
                        help='LSV IDs, separated by spaces, which should remain in the results. e.g. LSV_ID1 '
                             'GENE_ID2 ...')
filter_parser.add_argument('--lsv-ids-file', type=check_file,
                        help='Specify LSV IDs to retain in a file instead of on the command line. One LSV ID per line')
filter_parser.add_argument('--changing-threshold', type=float, default=0.2,
                        help='Threshold in delta-PSI quantification column. The default is "0.2".')
filter_parser.add_argument('--non-changing-threshold', type=float, default=0.05,
                        help='Threshold in delta-PSI quantification column. The default is "0.05".')
filter_parser.add_argument('--probability-changing-threshold', type=float, default=0.95,
                        help='The default is "0.95"')
filter_parser.add_argument('--probability-non-changing-threshold', type=float, default=0.95,
                        help='The default is "0.95"')
filter_parser.add_argument('--changing', action='store_true',
                         help='In general, find classifications for events which are changing (between multiple analysis). Requires at least one deltapsi or het voila file')
filter_parser.add_argument('--non-changing', action='store_true',
                         help='In general, find classifications for events which are not changing (between multiple analysis). Requires at least one deltapsi or het voila file')
required_filter_parser = filter_parser.add_argument_group('required named arguments')
required_filter_parser.add_argument('-d', '--directory', required=True, help="All filtered files will be dumped in"
                                                                          " this directory")


split_parser = argparse.ArgumentParser(add_help=False)
split_parser.add_argument('files', nargs='+', type=check_file,
                         help='List of files or directories which contains the splice graph and voila files.')
split_parser.add_argument('--copy-only', action='store_true',
                         help='The input files will not actually be split at all, they will just be duplicated to the '
                              'output directories as if they had been split. This may make the process much faster, '
                              'if you have the disk space to spare')
split_parser.add_argument('--overwrite', action='store_true',
                         help='If there are files inside of specified --directory, delete them and run splitter anyway')

required_split_parser = split_parser.add_argument_group('required named arguments')
required_split_parser.add_argument('-d', '--directory', required=True, help="Directories for each split will be created"
                                                                          " in this directory")
required_split_parser.add_argument('-n', '--num-divisions', required=True, type=int,
                         help='The number of separate directories to split the input into')


recombine_parser = argparse.ArgumentParser(add_help=False)
required_recombine_parser = recombine_parser.add_argument_group('required named arguments')
required_recombine_parser.add_argument('directories', nargs='+', help="Directory or directories"
                                                                " to search recursively for tsv files to combine")
required_recombine_parser.add_argument('-d', '--directory', required=True, help="Directory where recombined classifier "
                                                                                "output will be created.")


# subparsers
subparsers = parser.add_subparsers(help='')
subparsers.add_parser('tsv', parents=[tsv_parser, sys_parser, log_parser],
                      help='Generate tsv output for the supplied files.').set_defaults(func=Tsv)
subparsers.add_parser('view', parents=[view_parser, sys_parser, log_parser, webserver_parser],
                      help='Start service to view the visualization for the supplied files.').set_defaults(
    func=run_service)
subparsers.add_parser('modulize', parents=[classify_parser, sys_parser, log_parser],
                      help='Modulize splicing events and generate a breakdown of the modulization in '
                           'multiple TSV files.').set_defaults(func=Classify)
subparsers.add_parser('filter', parents=[filter_parser, sys_parser, log_parser],
                      help='Make truncated versions of the input files of a more manageable file size for easier'
                           ' collaboration.').set_defaults(
    func=Filter)
subparsers.add_parser('split', parents=[split_parser, sys_parser, log_parser],
                      help='Split classifier input dataset (splicegraph and voila files) into <N> equal parts'
                           ', for the purpose of running on a compute cluster').set_defaults(func=splitter)
subparsers.add_parser('recombine', parents=[recombine_parser, sys_parser, log_parser],
                      help='Used to combine output from a `voila split` run, after all initial '
                           'runs are complete').set_defaults(func=recombine)



if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

args = parser.parse_args()

if args.func == Classify:
    args.logger = args.directory + '/voila.log'


log = voila_log(filename=args.logger, silent=args.silent, debug=args.debug)

# dump all options on debug
for arg in vars(args):
    log.debug(f"Argument; {arg}: {getattr(args, arg)}")


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
