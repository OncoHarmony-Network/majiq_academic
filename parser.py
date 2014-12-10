import argparse

from pipelines import builder, calcpsi, deltapair

VERSION = "beta"


def new_subparser():
    return argparse.ArgumentParser(add_help=False)


def main():
    """
    Main MAJIQ parser with all flags and subcommands
    """
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



    buildparser = new_subparser()
    buildparser.add_argument('transcripts', action="store", help='read file in SAM format')
    buildparser.add_argument('-conf', default=None, help='Provide study configuration file with all '
                                                         'the execution information')
    buildparser.add_argument('-p', '--prefix', dest="prefix", type=str, default='', help='Output prefix string to '
                                                                                         'personalize partially the '
                                                                                         'output file.')
    buildparser.add_argument('--pcr', dest='pcr_filename', action="store", help='PCR bed file as gold_standard')
    buildparser.add_argument('--gff_output', dest='gff_output', default="lsvs.gff", action="store",
                             help='Filename where a gff with the lsv events will be generated')
    buildparser.add_argument('--minreads', default=2, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'the LSV "exist in the data". '
                             '[Default: %(default)s]')
    buildparser.add_argument('--minpos', default=2, type=int, help='Minimum number of start positions with at least 1 '
                                                                   'read in a LSV to consider that the LSV "exist in '
                                                                   'the data"')
<<<<<<< HEAD
    buildparser.add_argument('--only_rna', default=False, type=bool, action='store_true', help='Use only rna detected '
                                                                                               'junction in order to '
                                                                                               'detect LSV. If an exon '
                                                                                               'has only one junction '
                                                                                               'with coverage, it is '
                                                                                               'not going to be '
                                                                                               'detected as an LSV')
    buildparser.add_argument('--non_denovo', default=False, type=bool, action='store_true', help='Avoid denovo '
                                                                                                 'detection of '
                                                                                                 'junction, splicesites'
                                                                                                 ' and exons. This will'
                                                                                                 ' speedup the '
                                                                                                 'execution but reduce '
                                                                                                 'the number of LSVs '
                                                                                                 'detected')
=======
    buildparser.add_argument('--only_gather', action='store_true', dest='onlygather', default=False)
>>>>>>> master

    #flags shared by calcpsi and deltapair
    psianddelta = new_subparser()
    psianddelta.add_argument('--k', default=50, type=int,
                             help='Number of positions to sample per iteration. [Default: %(default)s]')
    psianddelta.add_argument('--m', default=100, type=int,
                             help='Number of bootstrapping samples. [Default: %(default)s]')
    psianddelta.add_argument('--minreads', default=10, type=int,
                             help='Minimum number of reads combining all positions in an event to be considered. '
                             '[Default: %(default)s]')
    psianddelta.add_argument('--minpos', default=3, type=int, help='Minimum number of start positions with at least 1 '
                                                                   'read for an event to be considered.')
    psianddelta.add_argument('--trimborder', default=5, type=int,
                             help='Trim the borders when sampling (keeping the ones with reads). '
                                  '[Default: %(default)s]')
    psianddelta.add_argument('--markstacks', default=0.0000001, type=float,
                             help='Mark stack positions. Expects a p-value. Use a negative value in order to '
                                  'disable it. [Default: %(default)s]')
    psianddelta.add_argument('--nogc', dest="gcnorm", action='store_false', default=True,
                             help='psianddelta GC content normalization [Default: GC content normalization activated]')
    psianddelta.add_argument('--nodiscardb', dest="discardb", action='store_false',  default=True,
                             help='Skip biscarding the b from the NB polynomial function, since we expect our fit '
                                  'to start from x=0, y=0')
    psianddelta.add_argument('--discardzeros', default=5, type=int, dest="discardzeros",
                             help='Discarding zeroes, up to a minimum of N positions per junction. [Default: 5]')

    #deltapair and deltagroup flags
    delta = new_subparser()
    delta.add_argument('-grp1', dest="files1", nargs='+', required=True)
    delta.add_argument('-grp2', dest="files2", nargs='+', required=True)
    delta.add_argument('--default_prior', action='store_true', default=False,
                       help="Use a default prior instead of computing it using the empirical data")
    delta.add_argument('--pairwise', default=False, action='store_true', help='')
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

    # delta.add_argument('--fixweights1', nargs='*', type=float,
    #                    help='Manually fix the weights for the replicas [Default: Automatic weight calculation]')
    # delta.add_argument('--fixweights2', nargs='*', type=float,
    #                    help='Manually fix the weights for the replicas [Default: Automatic weight calculation]')
    # delta.add_argument('--weightsL1', action='store_true', default=False, help='Use L1 instead of DKL in the weights '
    #                                                                            'algorithm')

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
                                              parents=[common, delta, psianddelta])
    parser_deltagroup.set_defaults(func=deltapair)
    args = parser.parse_args()
    args.func(args)



if __name__ == '__main__':
    main()






