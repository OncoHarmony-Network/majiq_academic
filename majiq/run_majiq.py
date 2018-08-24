import argparse
from majiq.src.build import build
from majiq.src.py_calc_psi import calcpsi as py_calcpsi
from majiq.src.py_deltapsi import deltapsi as py_deltapsi
from majiq.src.py_indpnt import calc_independent as py_calc_independent
from majiq.src.calc_psi import calcpsi
from majiq.src.deltapsi import deltapsi
from majiq.src.indpnt import calc_independent
from majiq.src.constants import *
from majiq.src.wght_pipeline import calc_weights as py_calc_weights

# from majiq.src.constants import all_stats
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


def new_subparser():
    return argparse.ArgumentParser(add_help=False)


def main():
    """
    Main MAJIQ parser with all flags and subcommands
    """
    #REMINDER parser.add_parser(..... parents='[bla, ble]')
    parser = argparse.ArgumentParser(description="MAJIQ is a suite of tools for the Splicing Events "
                                                 "and Alternative Splicing Quantification.")

    parser.add_argument('-v', action='version', version=VERSION)


    common = new_subparser()
    common.add_argument('-j', '--nproc', default=4, type=int, help='Number of processes to use')

    common.add_argument('-o', '--output', dest="outDir", required=True, help='Path to save the pickle output to.')
    common.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    common.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')

    common.add_argument('--debug', default=False, action='store_true',
                        help="Activate this flag for debugging purposes, activates logger and jumps some "
                             "processing steps. [Default: %(default)s]")
    common.add_argument('--mem-profile', default=False, action='store_true',
                        help="Print memory usage summary at the end of the execution. [Default: %(default)s]")

    common.add_argument('--min-experiments', default=-1, type=float, dest='min_exp',
                        help='Lower threshold for group filters. min_experiments is the minimum number of experiments '
                             'where the different filters check in order to pass an lsv or junction.')

    common.add_argument('--plotpath', default=None,
                        help='Path to save the plot to, if not provided will show on a matplotlib popup window')

    buildparser = new_subparser()
    buildparser.add_argument('transcripts', action="store", help='Annotation db ')
    buildparser.add_argument('-c', '--conf', default=None, required=True,
                             help='Provide study configuration file with all the execution information')

    buildparser.add_argument('--disable-ir', dest="ir", action='store_false', default=True,
                             help='Disables intron retention detection [Default: ir enabled]')

    buildparser.add_argument('--disable-denovo', dest="denovo", action='store_false', default=True,
                             help='Disables denovo detection of junction, splicesites and exons. This will speedup the '
                                  'execution but reduce the number of LSVs detected. [Default: denovo enabled]')

    buildparser.add_argument('--gff-output', dest='gff_output', default="lsvs.gff", action="store",
                             help='Filename where a gff with the lsv events will be generated. [Default: %(default)s]')

    buildparser.add_argument('--min-intronic-cov', default=1, type=float,
                             help='Minimum number of reads on average in intronic sites, only for intron retention.'
                                  'Default: %(default)s]')

    buildparser.add_argument('--min-denovo', default=2, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'denovo junction is real". [Default: %(default)s]')

    buildparser.add_argument('--minreads', default=3, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'the LSV "exist in the data". '
                             '[Default: %(default)s]')

    buildparser.add_argument('--minpos', default=2, type=int,
                             help='Minimum number of start positions with at least 1 read in a LSV to consider that '
                                  'the LSV "exist in the data". [Default: %(default)s]')

    # buildparser.add_argument('--only_rna', default=False, action='store_true',
    #                          help='Use only rna detected junction in order to detect LSV. If an exon has only one '
    #                               'junction with coverage, it is not going to be detected as an LSV. '
    #                               '[Default: %(default)s]')

    buildparser.add_argument('--use-db', default=False, action='store_true',
                             help='We provide a db.hdf5 file instead of a gff3 annotation db. This file is '
                                  'generated by a previous run of majiq [Default: %(default)s]')

    buildparser.add_argument('--markstacks', default=0.0000001, type=float, dest="pvalue_limit",
                             help='Mark stack positions. Expects a p-value. Use a negative value in order to '
                                  'disable it. [Default: %(default)s]')
    # buildparser.add_argument('--simplify', nargs='*')

    buildparser.add_argument('--k', default=50, type=int,
                             help='Number of positions to sample per iteration. [Default: %(default)s]')
    buildparser.add_argument('--m', default=30, type=int,
                             help='Number of bootstrapping samples. [Default: %(default)s]')


    sampling = new_subparser()

    sampling.add_argument('--minreads', default=10, type=int,
                          help='Minimum number of reads combining all positions in an event to be considered. '
                               '[Default: %(default)s]')
    sampling.add_argument('--minpos', default=3, type=int,
                          help='Minimum number of start positions with at least 1 read for an event to be considered.'
                               '[Default: %(default)s]')
    weights = new_subparser()
    weights.add_argument('--weights-alpha', type=float, default=15.,
                         help='Dispersion hyperparameter (Default: %(default)0.02f)')
    weights.add_argument('--weights-threshold', action=FRange01, default=0.75,
                         help='Threshold hyperparameter (Default: %(default)0.02f)')
    weights.add_argument('--weights-local', dest='local', type=float, default=0.,
                         help='Window for computation of local weights. If negative, uses a parametric approximation '
                              'instead. [Default: %(default)0.02f]')

    psi = new_subparser()
    psi.add_argument('files', nargs='+', help='The experiment files to analyze. You can include more than one '
                                              '(they will all be analyzed independently though) Glob syntax supported.')
    psi.add_argument('-n', '--name', required=True, help="The names that identify each of the experiments.")
    psi.add_argument('--weights', dest="weights", default='None',
                     help='Defines weights for each one of the replicas, for group1 and group2. The expected '
                          'value is --weights [Auto|None|<w1[,w2,..]>]\n'
                          '\t Auto will make majiq calculate the best weights, None will use uniform weights. '
                          'Select the weights manually requires specifying one weight for each replica or an '
                          'error will be triggered.')

    delta = new_subparser()
    delta.add_argument('-grp1', dest="files1", nargs='+', required=True)
    delta.add_argument('-grp2', dest="files2", nargs='+', required=True)
    delta.add_argument('--default-prior', action='store_true', default=False,
                       help="Use a default prior instead of computing it using the empirical data")
    delta.add_argument('-n', '--names', nargs='+', required=True,
                       help="The names that identify each of the experiments.")
    delta.add_argument('--binsize', default=0.025, type=int,
                       help='The bins for PSI values. With a --binsize of 0.025 (default), we have 40 bins')
    delta.add_argument('--prior-minreads', default=20, type=int,
                       help="Minimum number of reads combining all positions in a junction to be considered "
                            "(for the 'best set' calculation). [Default: %(default)s]")
    delta.add_argument('--prior-minnonzero', default=10, type=int,
                       help='Minimum number of positions for the best set.')
    delta.add_argument('--prior-iter', default=10, type=int, dest="iter",
                       help='Max number of iterations of the EM')

    delta.add_argument('--weights', dest="weights", nargs=2, default=['None', 'None'],
                       help='Defines weights for each one of the replicas, for group1 and group2. The expected '
                            'value is --weights [Auto|None|<w1[,w2,..]>] [Auto|None|<w\'1[,w\'2,..]>]\n'
                            '\t Auto will make majiq calculate the best weights, None will use uniform weights. '
                            'Select the weights manually requires specifying one weight for each replica or an '
                            'error will be triggered.')

    wght = new_subparser()
    wght.add_argument('files', nargs='+', help='The experiment files to analyze. You can include more than one '
                                               '(they will be analyzed independently though) Glob syntax supported.')
    wght.add_argument('-n', '--name', required=True, help="The names that identify each of the experiments. "
                                                          "[Default: %(default)s]")

    htrgen = new_subparser()
    htrgen.add_argument('-grp1', dest="files1", nargs='+', required=True)
    htrgen.add_argument('-grp2', dest="files2", nargs='+', required=True)
    htrgen.add_argument('-n', '--names', nargs='+', required=True,
                        help="The names that identify each of the experiments.")
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

    #calcpsi flags
    subparsers = parser.add_subparsers(help='')
    parser_preprocess = subparsers.add_parser('build', help='Preprocess SAM/BAM files as preparation for the rest of '
                                                            'the tools (psi, deltapsi)', parents=[common, buildparser])
    parser_preprocess.set_defaults(func=build)

    parser_calcpsi = subparsers.add_parser('py_psi', help="Calculate PSI values for N experiments, given a folder of "
                                                       "preprocessed events by 'majiq preprocess' or SAM/BAM files",
                                           parents=[common, psi, sampling, weights])
    parser_calcpsi.set_defaults(func=py_calcpsi)

    parser_deltagroup = subparsers.add_parser('py_deltapsi', help='Calculate Delta PSI values given a pair of experiments '
                                                               '(1 VS 1 conditions *with* replicas)',
                                              parents=[common, delta, sampling, weights])
    parser_deltagroup.set_defaults(func=py_deltapsi)

    parser_weights = subparsers.add_parser('weights', help='Calculate weights values given a group of experiment '
                                                           'replicas',
                                           parents=[common, sampling, weights, wght])
    parser_weights.set_defaults(func=py_calc_weights)

    parser_heterogen = subparsers.add_parser('py_heterogen', help='Calculate Delta PSI values given a pair of experiments '
                                                             'groups. This approach does not assume underlying PSI)',
                                             parents=[common, sampling, htrgen])
    parser_heterogen.set_defaults(func=py_calc_independent)


    parser_calcpsi = subparsers.add_parser('psi', help="Calculate PSI values for N experiments, given a folder of "
                                                       "preprocessed events by 'majiq preprocess' or SAM/BAM files",
                                           parents=[common, psi, sampling, weights])
    parser_calcpsi.set_defaults(func=calcpsi)

    parser_deltagroup = subparsers.add_parser('deltapsi', help='Calculate Delta PSI values given a pair of experiments '
                                                               '(1 VS 1 conditions *with* replicas)',
                                              parents=[common, delta, sampling, weights])
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






