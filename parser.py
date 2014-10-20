import argparse

from pipelines import builder, calcpsi, deltapair

VERSION = "alpha"


def new_subparser():
    return argparse.ArgumentParser(add_help=False)


def main():
    "Main MAJIQ parser with all flags and subcommands"
    #REMINDER parser.add_parser(..... parents='[bla, ble]')
    parser = argparse.ArgumentParser(description="MAJIQ is a suite of tools for the analysis of Alternative "
                                                 "Splicing Events and Alternative Splicing Quantification.")
    parser.add_argument('-v', action='version', version=VERSION)

    #common flags (first ones are required)
    common = new_subparser()
    common.add_argument('--nthreads', default=4, type=int, help='Number of threads')
    common.add_argument('--tmp', default="/tmp/", help='Path to save the temporary files. [Default: %(default)s]')
    common.add_argument('--output', required=True, help='Path to save the pickle output to.')
    common.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    common.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')
    common.add_argument('--plotpath', default=None,
                        help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    common.add_argument('--debug', type=int, default=0,
                        help="Activate this flag for debugging purposes, activates logger and jumps some "
                             "processing steps.")
    common.add_argument('--minreads', default=10, type=int,
                        help='Minimum number of reads combining all positions in an event to be considered. '
                             '[Default: %(default)s]')
    common.add_argument('--minnonzero', default=3, type=int, help='Minimum number of start positions with at least 1 '
                                                                  'read for an event to be considered.')


    buildparser = new_subparser()
    buildparser.add_argument('transcripts', action="store", help='read file in SAM format')
    buildparser.add_argument('-conf', default=None, help='Provide study configuration file with all '
                                                         'the execution information')
    buildparser.add_argument('-p', '--prefix', dest="prefix", type=str, default='', help='Output prefix string to '
                                                                                         'personalize partially the '
                                                                                         'output file.')
    buildparser.add_argument('--pcr', dest='pcr_filename', action="store", help='PCR bed file as gold_standard')
    buildparser.add_argument('--gff_output', dest='gff_output', action="store", help='Filename where a gff with the '
                                                                                     'lsv events will be generated')

    #flags shared by calcpsi and deltapair
    psianddelta = new_subparser()
    psianddelta.add_argument('--trim', default=0, type=int,
                             help='Trim the borders of the junctions because of poor mappability')
    psianddelta.add_argument('--k', default=50, type=int,
                             help='Number of positions to sample per iteration. [Default: %(default)s]')
    psianddelta.add_argument('--m', default=100, type=int,
                             help='Number of bootstrapping samples. [Default: %(default)s]')
    psianddelta.add_argument('--trimborder', default=5, type=int,
                             help='Trim the borders when sampling (keeping the ones with reads). '
                                  '[Default: %(default)s]')
    psianddelta.add_argument('--alpha', default=0.5, type=float,
                             help='Alpha hyperparameter for the dirichlet distribution. [Default: %(default)s]')
    psianddelta.add_argument('--markstacks', default=0.0000001, type=float,
                             help='Mark stack positions. Expects a p-value. Use a negative value in order to '
                                  'disable it. [Default: %(default)s]')
    psianddelta.add_argument('--nbdisp', default=0.1, type=float,
                             help='Dispersion for the fallback Negative Binomial function. Default: %(default)s]')
    psianddelta.add_argument('--nogc', dest="gcnorm", action='store_false', default=True,
                             help='psianddelta GC content normalization [Default: GC content normalization activated]')
    psianddelta.add_argument('--nodiscardb', dest="discardb", action='store_false',  default=True,
                             help='Skip biscarding the b from the NB polynomial function, since we expect our fit '
                                  'to start from x=0, y=0')
    psianddelta.add_argument('--discardzeros', default=5, type=int, dest="discardzeros",
                             help='Discarding zeroes, up to a minimum of N positions per junction. [Default: 5]')
    psianddelta.add_argument('--n', default=1, type=int,
                             help='Number of PSI samples per sample paired. [Default: %(default)s]')
    psianddelta.add_argument('--psiparam', default=False, action='store_true',
                             help='Instead of sampling, use a parametric form for the PSI calculation. '
                                  '[Default: %(default)s]')
    psianddelta.add_argument('--nz', default=0, type=int,
                             help='Method for number of non-zero position estimation.[0 - empirical, -1 - Binomial, '
                                  '>0 - Fized Value ]. Default: %(default)s]')

    #deltapair and deltagroup flags
    pairandgroup = new_subparser() 
    pairandgroup.add_argument('--names', nargs='+', required=True,
                              help="The names that identify each of the experiments. [Default: %(default)s]")
    pairandgroup.add_argument('--binsize', default=0.025, type=int,
                              help='The bins for PSI values. With a --binsize of 0.025 (default), we have 40 bins')
    pairandgroup.add_argument('--priorminreads', default=20, type=int,
                              help='Minimum number of reads combining all positions in a junction to be considered '
                                   '(for the "best set" calculation). [Default: %(default)s]')
    pairandgroup.add_argument('--priorminandreads', default=1, type=int,
                              help='Minimum number of reads combining all positions in a junction to be considered '
                                   '(for the "best set" calculation). [Default: %(default)s]')
    pairandgroup.add_argument('--priorminnonzero', default=10, type=int,
                              help='Minimum number of positions for the best set.')
    pairandgroup.add_argument('--iter', default=10, type=int,
                              help='Max number of iterations of the EM')
    pairandgroup.add_argument('--breakiter', default=0.01, type=float,
                              help='If the log likelihood increases less that this flag, do not do another EM step')
    pairandgroup.add_argument('--V', default=0.1, type=float,
                              help='Value of DeltaPSI used for initialization of the EM model [Default: %(default)s]')
    pairandgroup.add_argument('--synthprior', action='store_true', default=False,
                              help=' Generate the prior for DELTA PSI using our assumptions instead of the empirical '
                                   'data [Default: %(default)s]')
    pairandgroup.add_argument('--jefferiesprior', action='store_true', default=False,
                              help='Use only the jefferies prior, without including the  [Default: %(default)s]')
    pairandgroup.add_argument('--priorstd', default=0.15, type=float,
                              help="Standard deviation from the 0.5 PSI mean that the synthetic prior matrix has. "
                                   "Only works with --synthprior. [Default: %(default)s]")
    pairandgroup.add_argument('--prioruniform', default=3, type=float,
                              help="Uniform distribution to give a bit more of a chance to values out of the normal "
                                   "distribution. that the synthetic prior matrix has. Only works with --synthprior. "
                                   "[Default: %(default)s]")

    delta = new_subparser()
    delta.add_argument('-grp1', dest="files1", nargs='+')
    delta.add_argument('-grp2', dest="files2", nargs='+')
    delta.add_argument('--default_prior', action='store_true', default=False,
                       help="Use a default prior instead of computing it using the empirical data")
    delta.add_argument('--changinglimit')
    delta.add_argument('--changsetpercentile', type=float, default=90.,
                       help="Percentile of events that go into the 'best changing events' set")
    delta.add_argument('--fixweights1', nargs='*', type=float,
                       help='Manually fix the weights for the replicas [Default: Automatic weight calculation]')
    delta.add_argument('--fixweights2', nargs='*', type=float,
                       help='Manually fix the weights for the replicas [Default: Automatic weight calculation]')
    delta.add_argument('--weightsL1', action='store_true', default=False, help='Use L1 instead of DKL in the weights '
                                                                               'algorithm')
    delta.add_argument('--replicaweights', action='store_true', default=False,
                       help='Weight the experiments according to the events that change the most within replicas')
    delta.add_argument('--numbestchanging', default=200, type=int,
                       help="Number of events included in the best changing set (default %(default)s, should be "
                            "automatically calculated using FDR)")

    #calcpsi flags
    psi = new_subparser()
    psi.add_argument('files', nargs='+', help='The experiment files to analyze. You can include more than one '
                                              '(they will all be analyzed independently though) Glob syntax supported.')
    psi.add_argument('--name', required=True, help="The names that identify each of the experiments. "
                                                   "[Default: %(default)s]")


    subparsers = parser.add_subparsers(help='')
    parser_preprocess = subparsers.add_parser('build', help='Preprocess SAM/BAM files as preparation for the rest of '
                                                            'the tools (psi, deltapsi)', parents=[common, buildparser])
    parser_preprocess.set_defaults(func=builder)
    parser_calcpsi = subparsers.add_parser('psi', help="Calculate PSI values for N experiments, given a folder of "
                                                       "preprocessed events by 'majiq preprocess' or SAM/BAM files",
                                           parents=[common, psi, psianddelta])
    parser_calcpsi.set_defaults(func=calcpsi)
    parser_deltagroup = subparsers.add_parser('deltapsi', help='Calculate Delta PSI values given a pair of experiments '
                                                               '(1 VS 1 conditions *with* replicas)',
                                              parents=[common, delta, psianddelta, pairandgroup])
    parser_deltagroup.set_defaults(func=deltapair)
    args = parser.parse_args()
    args.func(args)



if __name__ == '__main__':
    main()






