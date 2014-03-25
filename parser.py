import argparse

from pipelines import preprocess, calcpsi, deltapair, deltagroup

VERSION = "alpha"


def new_subparser():
    return argparse.ArgumentParser(add_help=False)

def main():
    "Main MAJIQ parser with all flags and subcommands"
    #REMINDER parser.add_parser(..... parents='[bla, ble]')
    parser = argparse.ArgumentParser(description="MAJIQ is a suite of tools for the analysis of Alternative Splicing Events and Alternative Splicing Quantification.")
    parser.add_argument('-v', action='version', version=VERSION)

    #common flags (first ones are required)
    common = new_subparser()
    common.add_argument('--tmp', default="/tmp/", help='Path to save the temporary files. [Default: %(default)s]')
    common.add_argument('--output', required=True, help='Path to save the pickle output to.')
    common.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    common.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')
    common.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    common.add_argument('--debug', type=int, default=0, help="Activate this flag for debugging purposes, activates logger and jumps some processing steps.")

    #flags shared by calcpsi and deltapair
    psianddelta = new_subparser()
    psianddelta.add_argument('--trim', default=0, type=int, help='Trim the borders of the junctions because of poor mappability')
    psianddelta.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration. [Default: %(default)s]')
    psianddelta.add_argument('--m', default=100, type=int, help='Number of bootstrapping samples. [Default: %(default)s]')  
    psianddelta.add_argument('--trimborder', default=5, type=int, help='Trim the borders when sampling (keeping the ones with reads). [Default: %(default)s]')
    psianddelta.add_argument('--alpha', default=0.5, type=float, help='Alpha hyperparameter for the dirichlet distribution. [Default: %(default)s]') 
    psianddelta.add_argument('--markstacks', default=0.001, type=float, help='Mark stack positions. Expects a p-value. Use a negative value in order to disable it. [Default: %(default)s]') 
    psianddelta.add_argument('--nbdisp', default=0.1, type=int, help='Dispersion for the fallback Negative Binomial function. [Default: %(default)s]')
    psianddelta.add_argument('--nogc', dest="gcnorm", action='store_false', default=True, help='psianddelta GC content normalization [Default: GC content normalization activated]')
    psianddelta.add_argument('--nodiscardb', dest="discardb", action='store_false',  default=True, help='Skip biscarding the b from the NB polynomial function, since we expect our fit to start from x=0, y=0')
    psianddelta.add_argument('--nodiscardzeros', action='store_false', default=True, dest="discardzeros", help='Skip discarding zeroes')    
    psianddelta.add_argument('--n', default=1, type=int, help='Number of PSI samples per sample paired. [Default: %(default)s]') 
    psianddelta.add_argument('--psiparam', default=False, action='store_true', help='Instead of sampling, use a parametric form for the PSI calculation. [Default: %(default)s]')

    #psianddelta.add_argument('--ONLYSTACKS', action='store_true', help="Debug flag that should dissapear. Used to test if stacks are worth masking.")
    #psianddelta.add_argument('--usetensor', action='store_true')

    #deltapair flags
    delta = new_subparser()
    delta.add_argument('file1')  
    delta.add_argument('file2') 

    #deltapair and deltagroup flags
    pairandgroup = new_subparser() 
    pairandgroup.add_argument('--names', nargs='+', required=True, help="The names that identify each of the experiments. [Default: %(default)s]")
    pairandgroup.add_argument('--binsize', default=0.025, type=int, help='The bins for PSI values. With a --binsize of 0.025 (default), we have 40 bins')   
    pairandgroup.add_argument('--minreads', default=20, type=int, help='Minimum number of reads combining all positions in a junction to be considered (for the "best set" calculation). [Default: %(default)s]') 
    pairandgroup.add_argument('--minandreads', default=1, type=int, help='Minimum number of reads combining all positions in a junction to be considered (for the "best set" calculation). [Default: %(default)s]') 
    pairandgroup.add_argument('--minnonzero', default=10, type=int, help='Minimum number of positions for the best set.')
    pairandgroup.add_argument('--iter', default=10, type=int, help='Max number of iterations of the EM')
    pairandgroup.add_argument('--breakiter', default=0.01, type=float, help='If the log likelihood increases less that this flag, do not do another EM step')
    pairandgroup.add_argument('--V', default=0.1, type=float, help='Value of DeltaPSI used for initialization of the EM model [Default: %(default)s]')
    pairandgroup.add_argument('--synthprior', action='store_true', default=False, help=' Generate the prior for DELTA PSI using our assumptions instead of the empirical data [Default: %(default)s]')
    pairandgroup.add_argument('--jefferiesprior', action='store_true', default=False, help='Use only the jefferies prior, without including the  [Default: %(default)s]')
    pairandgroup.add_argument('--priorstd', default=0.15, type=float, help="Standard deviation from the 0.5 PSI mean that the synthetic prior matrix has. Only works with --synthprior. [Default: %(default)s]")
    pairandgroup.add_argument('--prioruniform', default=3, type=float, help="Uniform distribution to give a bit more of a chance to values out of the normal distribution. that the synthetic prior matrix has. Only works with --synthprior. [Default: %(default)s]")


    group = new_subparser()
    group.add_argument('--1', dest="files1", nargs='+')  
    group.add_argument('--2', dest="files2", nargs='+')
    group.add_argument('--changinglimit')
    group.add_argument('--changsetpercentile', type=float, default=90., help="Percentile of events that go into the 'best changing events' set")
    group.add_argument('--fixweights1', nargs='*', type=float, help='Manually fix the weights for the replicas [Default: Automatic weight calculation]')
    group.add_argument('--fixweights2', nargs='*', type=float, help='Manually fix the weights for the replicas [Default: Automatic weight calculation]')
    group.add_argument('--weightsL1', action='store_true', default=False, help='Use L1 instead of DKL in the weights algorithm')    
    group.add_argument('--replicaweights', action='store_true', default=False, help='Weight the experiments according to the events that change the most within replicas')
    group.add_argument('--numbestchanging', default=200, type=int, help="Number of events included in the best changing set (default %(default)s, should be automatically calculated using FDR)")

    #calcpsi flags
    psi = new_subparser()
    psi.add_argument('files', nargs='+', help='The experiment files to analyze. You can include more than one (they will all be analyzed independently though) Glob syntax supported.')
    psi.add_argument('--minreads', default=0, type=int, help='Minimum number of reads combining all positions in an event to be considered. [Default: %(default)s]') 
    psi.add_argument('--minnonzero', default=0, type=int, help='Minimum number of start positions with at least 1 read for an event to be considered.')

    subparsers = parser.add_subparsers(help='')
    parser_preprocess = subparsers.add_parser('preprocess', help='Preprocess SAM/BAM files as preparation for the rest of the tools (calcpsi, deltapair, deltagroup)', parents=[common])
    parser_preprocess.set_defaults(func=preprocess)
    parser_calcpsi = subparsers.add_parser('calcpsi', help="Calculate PSI values for N experiments, given a folder of preprocessed events by 'majiq preprocess' or SAM/BAM files (This last not implemented yet)", parents=[common, psi, psianddelta])
    parser_calcpsi.set_defaults(func=calcpsi)
    parser_deltapair = subparsers.add_parser('deltapair', help='Calculate Delta PSI values given a pair of experiments (1 VS 1 conditions without replicas)', parents=[common, delta, psianddelta, pairandgroup])
    parser_deltapair.set_defaults(func=deltapair)
    parser_deltagroup = subparsers.add_parser('deltagroup', help='Calculate Delta PSI values given a pair of experiments (1 VS 1 conditions *with* replicas)', parents=[common, group, psianddelta, pairandgroup])
    parser_deltagroup.set_defaults(func=deltagroup)
    args = parser.parse_args()
    args.func(args)



if __name__ == '__main__':
    main()






