import argparse
from majiq.src.build import build
from majiq.src.calc_psi import calcpsi
from majiq.src.deltapsi import deltapsi
from majiq.src.indpnt import calc_independent
from majiq.src.constants import *
import sys


class FRange01(argparse.Action):
    def __call__(self, parser, namespace, values, option_string = None):
        values = float(values)
        if values < 0 or values > 1:
            raise ValueError('must be in range [0, 1]')
        setattr(namespace, self.dest, values)


def get_vwindow(param):
    param = float(param)
    assert 0.0 <= param <= 1.0
    return param


def get_bins(param):
    return np.linspace(0, 1, num=int(param) + 1)


def get_minsamps(param):
    param = int(param)
    assert param >= 2
    return param

def check_positive(value):
    ivalue = float(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive float value" % value)
    return ivalue



def new_subparser():
    return argparse.ArgumentParser(add_help=False)


def main():
    """
    Main MAJIQ parser with all flags and subcommands
    """
    #REMINDER parser.add_parser(..... parents='[bla, ble]')
    parser = argparse.ArgumentParser(description="MAJIQ is a suite of tools for the Splicing Events "
                                                 "and Alternative Splicing Quantification.")

    parser.add_argument('-v', action='version', version="%s-%s" % (VERSION, get_git_version()))

    common = new_subparser()
    common.add_argument('-j', '--nproc', default=4, type=int, help='Number of processes to use. [Default: %(default)s]')
    common.add_argument('-o', '--output', dest="outDir", required=True, help='Path to save the pickle output to.')

    common.add_argument('--logger', default=None, help='Path for the logger. [Default is output directory]')
    common.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')

    common.add_argument('--debug', default=False, action='store_true',
                        help="Activate this flag for debugging purposes, activates logger and jumps some "
                             "processing steps. [Default: %(default)s]")

    common.add_argument('--mem-profile', default=False, action='store_true',
                        help="Print memory usage summary at the end of the execution. [Default: %(default)s]")

    common.add_argument('--min-experiments', default=0.5, type=check_positive, dest='min_exp',
                        help='Lower threshold for group filters. min_experiments set the minimum number of experiments '
                             'where the different filters check in order to pass an lsv or junction.\n'
                             '\t + <  1 the value is the fraction of the experiments in the group\n'
                             '\t + >= 1 the value is the actual number of experiments. If the number is set to a '
                             'greater number than the size of the group, we use the size instead.\n'
                             '[Default: %(default)s]]')

    common.add_argument('--plotpath', default=None,
                        help='Path to save the plot to, if not provided will show on a matplotlib popup window. '
                             '[Default: %(default)s]')

    buildparser = new_subparser()
    buildparser.add_argument('transcripts', action="store", help='Annotation db ')
    buildparser.add_argument('-c', '--conf', default=None, required=True,
                             help='Provide study configuration file with all the execution information')

    buildparser.add_argument('--disable-ir', dest="ir", action='store_false', default=True,
                             help='Disables intron retention detection [Default: ir enabled]')

    buildparser.add_argument('--disable-denovo', dest="denovo", action='store_false', default=True,
                             help='Disables denovo detection of junction, splicesites and exons. This will speedup the '
                                  'execution but reduce the number of LSVs detected. [Default: denovo enabled]')

    buildparser.add_argument('--min-intronic-cov', default=0.01, type=float,
                             help='Minimum number of reads on average in intronic sites, only for intron retention.'
                                  'Default: %(default)s]')

    buildparser.add_argument('--min-denovo', default=10, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'denovo junction is real". [Default: %(default)s]')

    buildparser.add_argument('--minreads', default=3, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'the LSV "exist in the data". '
                             '[Default: %(default)s]')

    buildparser.add_argument('--minpos', default=2, type=int,
                             help='Minimum number of start positions with at least 1 read in a LSV to consider that '
                                  'the LSV "exist in the data". [Default: %(default)s]')

    buildparser.add_argument('--markstacks', default=0.0000001, type=float, dest="pvalue_limit",
                             help='Mark stack positions. Expects a p-value. Use a negative value in order to '
                                  'disable it. [Default: %(default)s]')

    buildparser.add_argument('--k', default=50, type=int,
                             help='Number of positions to sample per iteration. [Default: %(default)s]')
    buildparser.add_argument('--m', default=30, type=int,
                             help='Number of bootstrapping samples. [Default: %(default)s]')
    buildparser.add_argument('--irnbins', default=0.5, type=float, help='This values defines the number of bins with '
                                                                        'some coverage that an intron needs to pass '
                                                                        'to be accepted as real [Default: %(default)s]')

    buildparser.add_argument('--simplify_denovo', dest="simpl_denovo", default=10, type=int,
                             help='Minimum number of reads threshold combining all positions of an denovo junction to '
                                  'consider if it will be simplified, even knowing it is real. Simplified junctions are'
                                  ' discarded from any lsv. [Default: %(default)s]')

    buildparser.add_argument('--simplify_annotated', dest="simpl_db", default=10, type=int,
                             help='Minimum number of reads threshold combining all positions of an annotated junction to '
                                  'consider if it will be simplified, even knowing it is real. Simplified junctions are'
                                  ' discarded from any lsv. [Default: %(default)s]')

    buildparser.add_argument('--simplify_ir', dest="simpl_ir", default=10, type=int,
                             help='Minimum number of reads threshold combining all positions of an ir to '
                                  'consider if it will be simplified, even knowing it is real. Simplified junctions are'
                                  ' discarded from any lsv. [Default: %(default)s]')

    buildparser.add_argument('--simplify_psi', dest="simpl_psi", default=0.01, type=float,
                             help='Minimum fraction of the usage of any junction in a LSV to consider that junction is'
                                  ' real. [Default: %(default)s]')



    sampling = new_subparser()

    sampling.add_argument('--minreads', default=10, type=int,
                          help='Minimum number of reads combining all positions in an event to be considered. '
                               '[Default: %(default)s]')
    sampling.add_argument('--minpos', default=3, type=int,
                          help='Minimum number of start positions with at least 1 read for an event to be considered.'
                               '[Default: %(default)s]')

    psi = new_subparser()
    psi.add_argument('files', nargs='+', help='The experiment files to analyze. You can include more than one '
                                              '(they will all be analyzed independently though) Glob syntax supported.')
    psi.add_argument('-n', '--name', required=True, help="The names that identify each of the experiments.")

    psi.add_argument('--output-type', choices=['voila', 'tsv', 'all'], default='all',
                     help='Defines the type of output file to be generated, voila file to be used in voila, '
                          'tsv with the lsv information or both. [Default: %(default)s]')

    delta = new_subparser()
    delta.add_argument('-grp1', dest="files1", nargs='+', required=True)
    delta.add_argument('-grp2', dest="files2", nargs='+', required=True)
    delta.add_argument('--default-prior', action='store_true', default=False,
                       help="Use a default prior instead of computing it using the empirical data. "
                            " [Default: default prior disabled]")
    delta.add_argument('-n', '--names', nargs=2, required=True,
                       help="The names that identify each of the experiments.")
    delta.add_argument('--binsize', default=0.025, type=int,
                       help='The bins for PSI values. With a --binsize of 0.025 (default), we have 40 bins. '
                            '[Default: %(default)s]')
    delta.add_argument('--prior-minreads', default=20, type=int,
                       help="Minimum number of reads combining all positions in a junction to be considered "
                            "(for the 'best set' calculation). [Default: %(default)s]")
    delta.add_argument('--prior-minnonzero', default=10, type=int,
                       help='Minimum number of positions for the best set. [Default: %(default)s]')
    delta.add_argument('--prior-iter', default=1, type=int, dest="iter",
                       help='Max number of iterations of the EM. [Default: %(default)s]')
    delta.add_argument('--output-type', choices=['voila', 'tsv', 'all'], default='all',
                       help='Defines the type of output file to be generated, voila file to be used in voila, '
                            'tsv with the lsv information or both. [Default: %(default)s]')

    htrgen = new_subparser()
    htrgen.add_argument('-grp1', dest="files1", nargs='+', required=True)
    htrgen.add_argument('-grp2', dest="files2", nargs='+', required=True)
    htrgen.add_argument('-n', '--names', nargs='+', required=True,
                        help="The names that identify each of the experiments.")
    htrgen.add_argument('--keep-tmpfiles', action='store_true', default=False, dest='keep_tmpfiles',
                        help='When this argument is specified, majiq heterogen will not remove the psi files that '
                             'are temporary generated during the execution [Default: %(default)d]')
    htrgen.add_argument('--nsamples', type=int, default=100, dest="psi_samples",
                        help='Number of PSI samples to take per LSV junction. If equal to 1, use expected value only. '
                             '[Default: %(default)d]')
    htrgen.add_argument('--vwindow', type=get_vwindow, default=0.0,
                        help='Width of sample rejection window. If equal to 0, do not reject samples. '
                             '[Default: %(default)0.02f.]')
    htrgen.add_argument('--bins', type=get_bins, default=get_bins(40),
                        help='Fixed-width binning resolution of PSI distributions. [Default: 40')
    htrgen.add_argument('--stats', nargs='+', default=['ALL'],
                        help='Test statistics to run. [Default: %(default)s]')
    htrgen.add_argument('--minsamps', type=get_minsamps, default=2,
                        help='Minimum number of samples that need to be present for an LSV junction in order to '
                             'perform each test statistic. [Default: %(default)d]')
    htrgen.add_argument('--test_percentile', type=float, default=0.95,
                        help='For each one of the statistical tests, we combine all pvalue per psi sample by '
                             'percentile calculation. This argument allows the user define with percentile they '
                             'want to use [Default: %(default)d]')





    #calcpsi flags
    subparsers = parser.add_subparsers(help='')
    parser_preprocess = subparsers.add_parser('build', help='Preprocess SAM/BAM files as preparation for the rest of '
                                                            'the tools (psi, deltapsi)', parents=[common, buildparser])
    parser_preprocess.set_defaults(func=build)

    parser_calcpsi = subparsers.add_parser('psi', help="Calculate PSI values for N experiments, given a folder of "
                                                       "preprocessed events by 'majiq preprocess' or SAM/BAM files",
                                           parents=[common, psi, sampling])
    parser_calcpsi.set_defaults(func=calcpsi)

    parser_deltagroup = subparsers.add_parser('deltapsi', help='Calculate Delta PSI values given a pair of experiments '
                                                               '(1 VS 1 conditions *with* replicas)',
                                              parents=[common, delta, sampling])
    parser_deltagroup.set_defaults(func=deltapsi)

    parser_heterogen = subparsers.add_parser('heterogen', help='Calculate Delta PSI values given a pair of experiments '
                                                             'groups. This approach does not assume underlying PSI)',
                                             parents=[common, sampling, htrgen])
    parser_heterogen.set_defaults(func=calc_independent)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()






