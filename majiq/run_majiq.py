import argparse
from majiq.src.build import build
from majiq.src.calc_psi import calcpsi
from majiq.src.deltapsi import deltapsi
from majiq.src.constants import *
from majiq.src.wght_pipeline import calc_weights


class FRange01(argparse.Action):
    def __call__(self, parser, namespace, values, option_string = None):
        values = float(values)
        if values < 0 or values > 1:
            raise ValueError('must be in range [0, 1]')
        setattr(namespace, self.dest, values)


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

    #common flags (first ones are required)
    common = new_subparser()
    common.add_argument('--nthreads', default=4, type=int, help='Number of threads')
    common.add_argument('--tmp', default="/tmp/", help='Path to save the temporary files. [Default: %(default)s]')
    common.add_argument('--output', dest="outDir", required=True, help='Path to save the pickle output to.')
    common.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    common.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')
    common.add_argument('--plotpath', default=None,
                        help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    common.add_argument('--debug', type=int, default=0,
                        help="Activate this flag for debugging purposes, activates logger and jumps some "
                             "processing steps.")
    common.add_argument('--min_experiments', default=-1, type=float, dest='min_exp')

    buildparser = new_subparser()
    buildparser.add_argument('transcripts', action="store", help='read file in SAM format')
    buildparser.add_argument('-conf', default=None, required=True, help='Provide study configuration file with all '
                                                                        'the execution information')
    buildparser.add_argument('--nogc', dest="gcnorm", action='store_false', default=True,
                             help='psianddelta GC content normalization [Default: GC content normalization activated]')
    buildparser.add_argument('--pcr', dest='pcr_filename', action="store", help='PCR bed file as gold_standard')
    buildparser.add_argument('--gff_output', dest='gff_output', default="lsvs.gff", action="store",
                             help='Filename where a gff with the lsv events will be generated')
    buildparser.add_argument('--min_denovo', default=2, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'denovo junction is real". '
                             '[Default: %(default)s]')
    buildparser.add_argument('--minreads', default=3, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'the LSV "exist in the data". '
                             '[Default: %(default)s]')
    buildparser.add_argument('--min_intronic_cov', default=1.5, type=float,
                             help='Minimum number of reads on average in intronic sites, only for intron retention.'
                                  'Default: %(default)s]')
    buildparser.add_argument('--minpos', default=2, type=int, help='Minimum number of start positions with at least 1 '
                                                                   'read in a LSV to consider that the LSV "exist in '
                                                                   'the data"')

    buildparser.add_argument('--only_rna', default=False, action='store_true', help='Use only rna detected junction in '
                                                                                    'order to detect LSV. If an exon '
                                                                                    'has only one junction with '
                                                                                    'coverage, it is not going to be '
                                                                                    'detected as an LSV')
    buildparser.add_argument('--non_denovo', default=False, action='store_true', help='Avoid denovo detection of '
                                                                                      'junction, splicesites and exons.'
                                                                                      ' This will speedup the execution'
                                                                                      ' but reduce the number of LSVs '
                                                                                      'detected')
    buildparser.add_argument('--only_gather', action='store_true', dest='onlygather', default=False)
    buildparser.add_argument('--permissive_ir', action='store_true', dest='permissive_ir', default=False)
    buildparser.add_argument('--markstacks', default=0.0000001, type=float,
                             help='Mark stack positions. Expects a p-value. Use a negative value in order to '
                                  'disable it. [Default: %(default)s]')
    buildparser.add_argument('--simplify', nargs='*')

    buildparser.add_argument('--prebam', default=True,  action='store_false')

    sampling = new_subparser()
    sampling.add_argument('--k', default=50, type=int,
                             help='Number of positions to sample per iteration. [Default: %(default)s]')
    sampling.add_argument('--m', default=100, type=int,
                             help='Number of bootstrapping samples. [Default: %(default)s]')
    sampling.add_argument('--minreads', default=10, type=int,
                             help='Minimum number of reads combining all positions in an event to be considered. '
                                  '[Default: %(default)s]')
    sampling.add_argument('--minpos', default=3, type=int, help='Minimum number of start positions with at least 1 '
                                                                   'read for an event to be considered.'
                                                                   '[Default: %(default)s]')
    sampling.add_argument('--trimborder', default=5, type=int,
                             help='Trim the borders when sampling (keeping the ones with reads). '
                                  '[Default: %(default)s]')
    sampling.add_argument('--nodiscardb', dest="discardb", action='store_false', default=True,
                             help='Skip biscarding the b from the NB polynomial function, since we expect our fit '
                                  'to start from x=0, y=0')
    sampling.add_argument('--discardzeros', default=5, type=int, dest="discardzeros",
                             help='Discarding zeroes, up to a minimum of N positions per junction. [Default: 5]')


    #flags shared by calcpsi and deltapair
    weights = new_subparser()
    weights.add_argument('--weights_alpha', type=float, default=15.,
                             help='Dispersion hyperparameter (Default: %(default)0.02f)')
    weights.add_argument('--weights_threshold', action=FRange01, default=0.75,
                             help='Threshold hyperparameter (Default: %(default)0.02f)')
    weights.add_argument('--weights_local', dest='local', type=float, default=0.,
                            help='Window for computation of local weights. '
                                 'If negative, uses a parametric '
                                 'approximation instead. ('
                                 'Default: %(default)0.02f)')

    psi = new_subparser()
    psi.add_argument('files', nargs='+', help='The experiment files to analyze. You can include more than one '
                                              '(they will all be analyzed independently though) Glob syntax supported.')
    psi.add_argument('--name', required=True, help="The names that identify each of the experiments. "
                                                   "[Default: %(default)s]")
    psi.add_argument('--only_bootstrap', action='store_true', dest='only_boots', default=False)
    psi.add_argument('--weights', dest="weights", default='None',
                     help='Defines weights for each one of the replicas, for group1 and group2. The expected '
                          'value is --weights [Auto|None|<w1[,w2,..]>]\n'
                          '\t Auto will make majiq calculate the best weights, None will use uniform weights. '
                          'Select the weights manually requires specifying one weight for each replica or an '
                          'error will be triggered.')



    #deltapair and deltagroup flags
    delta = new_subparser()
    delta.add_argument('-grp1', dest="files1", nargs='+', required=True)
    delta.add_argument('-grp2', dest="files2", nargs='+', required=True)
    delta.add_argument('--default_prior', action='store_true', default=False,
                       help="Use a default prior instead of computing it using the empirical data")
    delta.add_argument('--names', nargs='+', required=True,
                       help="The names that identify each of the experiments. [Default: %(default)s]")
    delta.add_argument('--binsize', default=0.025, type=int,
                       help='The bins for PSI values. With a --binsize of 0.025 (default), we have 40 bins')
    delta.add_argument('--priorminreads', default=20, type=int,
                       help="Minimum number of reads combining all positions in a junction to be considered "
                            "(for the 'best set' calculation). [Default: %(default)s]")
    delta.add_argument('--priorminnonzero', default=10, type=int,
                       help='Minimum number of positions for the best set.')
    delta.add_argument('--iter', default=10, type=int,
                       help='Max number of iterations of the EM')
    delta.add_argument('--breakiter', default=0.01, type=float,
                       help='If the log likelihood increases less that this flag, do not do another EM step')
    delta.add_argument('--prioruniform', default=3, type=float,
                       help="Uniform distribution to give a bit more of a chance to values out of the normal "
                            "distribution. that the synthetic prior matrix has. Only works with --synthprior. "
                            "[Default: %(default)s]")
    delta.add_argument('--extra_bootstrap', action='store_true', dest='export_boots', default=False)
    delta.add_argument('--weights', dest="weights", nargs=2, default=['None', 'None'],
                       help='Defines weights for each one of the replicas, for group1 and group2. The expected '
                            'value is --weights [Auto|None|<w1[,w2,..]>] [Auto|None|<w\'1[,w\'2,..]>]\n'
                            '\t Auto will make majiq calculate the best weights, None will use uniform weights. '
                            'Select the weights manually requires specifying one weight for each replica or an '
                            'error will be triggered.')

    wght = new_subparser()
    wght.add_argument('files', nargs='+', help='The experiment files to analyze. You can include more than one '
                                              '(they will all be analyzed independently though) Glob syntax supported.')
    wght.add_argument('--name', required=True, help="The names that identify each of the experiments. "
                                                   "[Default: %(default)s]")


    # mdelta = new_subparser()
    # mdelta.add_argument('-pairs_file', dest="deltapairs", required=True)
    # mdelta.add_argument('--default_prior', action='store_true', default=False,
    #                    help="Use a default prior instead of computing it using the empirical data")
    # mdelta.add_argument('--binsize', default=0.025, type=int,
    #                    help='The bins for PSI values. With a --binsize of 0.025 (default), we have 40 bins')
    # mdelta.add_argument('--priorminreads', default=20, type=int,
    #                    help="Minimum number of reads combining all positions in a junction to be considered "
    #                         "(for the 'best set' calculation). [Default: %(default)s]")
    # mdelta.add_argument('--priorminnonzero', default=10, type=int,
    #                    help='Minimum number of positions for the best set.')
    # mdelta.add_argument('--iter', default=10, type=int,
    #                    help='Max number of iterations of the EM')
    # mdelta.add_argument('--breakiter', default=0.01, type=float,
    #                    help='If the log likelihood increases less that this flag, do not do another EM step')
    # mdelta.add_argument('--prioruniform', default=3, type=float,
    #                    help="Uniform distribution to give a bit more of a chance to values out of the normal "
    #                         "distribution. that the synthetic prior matrix has. Only works with --synthprior. "
    #                         "[Default: %(default)s]")

    #calcpsi flags
    subparsers = parser.add_subparsers(help='')
    parser_preprocess = subparsers.add_parser('build', help='Preprocess SAM/BAM files as preparation for the rest of '
                                                            'the tools (psi, deltapsi)', parents=[common, buildparser])
    parser_preprocess.set_defaults(func=build)

    parser_calcpsi = subparsers.add_parser('psi', help="Calculate PSI values for N experiments, given a folder of "
                                                       "preprocessed events by 'majiq preprocess' or SAM/BAM files",
                                           parents=[common, psi, sampling, weights])
    parser_calcpsi.set_defaults(func=calcpsi)

    parser_deltagroup = subparsers.add_parser('deltapsi', help='Calculate Delta PSI values given a pair of experiments '
                                                               '(1 VS 1 conditions *with* replicas)',
                                              parents=[common, delta, sampling, weights])
    parser_deltagroup.set_defaults(func=deltapsi)

    parser_weights = subparsers.add_parser('weights', help='Calculate weights values given a group of experiment '
                                                           'replicas',
                                              parents=[common, sampling, weights, wght])
    parser_weights.set_defaults(func=calc_weights)



    # parser_multidelta = subparsers.add_parser('multi_delta', help='Calculate Delta PSI values given a pair of experiments '
    #                                                            '(1 VS 1 conditions *with* replicas)',
    #                                           parents=[common, mdelta, psianddelta])
    # parser_multidelta.set_defaults(func=multi_dpsi)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()






