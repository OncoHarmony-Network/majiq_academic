#!/usr/bin/env python
"""
run_majiq.py

Entrypoint into different MAJIQ tools. Parse options from command-line and
launch the associated MAJIQ pipeline.

Authors: Jorge Vaquero-Garcia, Juan Gonzalez-Vallinas, Joseph K Aicher
"""

import argparse
from rna_majiq.src.build import build
from rna_majiq.src.calc_psi import calcpsi
from rna_majiq.src.deltapsi import deltapsi
from rna_majiq.src.indpnt import calc_independent
from rna_majiq.majiq_args import (
    StoreRequiredUniqueActionFactory,
    check_nonnegative_factory,
    check_characters_factory,
)
import rna_majiq.src.constants as constants
import sys
from rna_voila.api.licensing import check_license


StoreGroupNames = StoreRequiredUniqueActionFactory()
StoreExperimentPaths = StoreRequiredUniqueActionFactory()

check_positive = check_nonnegative_factory(float, True)
check_nonnegative = check_nonnegative_factory(float, False)
check_positive_int = check_nonnegative_factory(int, True)
check_nonnegative_int = check_nonnegative_factory(int, False)

check_group_chars = check_characters_factory(
    constants.ALLOWED_GROUP_NAME_CHARS, "alphanumeric or underscore characters"
)


def new_subparser():
    """ Just get an empty subparser
    """
    return argparse.ArgumentParser(add_help=False)


def main():
    """
    Main MAJIQ parser with all flags and subcommands
    """
    # REMINDER parser.add_parser(..... parents='[bla, ble]')
    parser = argparse.ArgumentParser(
        description="MAJIQ is a suite of tools to detect and quantify local"
        " splicing variations (LSVs) from RNA-seq data."
    )

    parser.add_argument('--license', required=False)

    parser.add_argument(
        "-v",
        action="version",
        version=f"{constants.VERSION}",
    )

    common = new_subparser()
    common.add_argument(
        "-j",
        "--nproc",
        default=4,
        type=check_positive_int,
        help="Number of threads to use. [Default: %(default)s]",
    )
    common.add_argument(
        "-o",
        "--output",
        dest="outDir",
        required=True,
        help="Path for output directory to which output files will be saved.",
    )

    common.add_argument(
        "--logger",
        default=None,
        help="Path for the logger. [Default is output directory]",
    )
    common.add_argument(
        "--silent", action="store_true", default=False, help="Silence the logger."
    )

    common.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="This flag is used for debugging purposes. It activates more"
        " verbose logging and skips some processing steps. [Default: %(default)s]",
    )

    common.add_argument(
        "--mem-profile",
        default=False,
        action="store_true",
        help="Print memory usage summary at the end of program execution."
        " [Default: %(default)s]",
    )

    common.add_argument(
        "--min-experiments",
        default=0.5,
        type=check_positive,
        dest="min_exp",
        help="Threshold for group filters. This specifies the fraction"
        " (value < 1) or absolute number (value >= 1) of experiments passing"
        " per-experiment filters (i.e. minreads, minpos, etc.) that must pass"
        " individually in order to pass an LSV or junction. If greater than"
        " the total number of experiments in a group, requires all experiments"
        " to pass individually. [Default: %(default)s]",
    )

    buildparser = new_subparser()
    buildparser_required = buildparser.add_argument_group("Required arguments")
    buildparser_required.add_argument(
        "transcripts", action="store", help="Annotation database in GFF3 format",
    )
    buildparser_required.add_argument(
        "-c",
        "--conf",
        default=None,
        required=True,
        help="Path to study configuration file specifying paths to BAM/SJ"
        " files, samples for each group, and sample-specific information."
        " Format follows rules for INI files, with required/optional fields"
        " used by MAJIQ described in the MAJIQ quick start:"
        " <https://biociphers.bitbucket.io/majiq/quick.html>",
    )

    # per-experiment filters for all junctions
    buildparser_junctions = buildparser.add_argument_group("Junction filters")
    buildparser_junctions.add_argument(
        "--minreads",
        default=3,
        type=check_positive_int,
        help="Threshold on the minimum total number of reads for any junction"
        " to meet per-experiment filters for the LSVs it is a part of. When"
        " the minimum numbers of reads and positions (--minpos) are both met"
        " in enough experiments for a group in any junction that is part of an"
        " LSV, the LSV is considered admissible and saved in output MAJIQ"
        " files for potential downstream quantification. [Default: %(default)s]",
    )
    buildparser_junctions.add_argument(
        "--minpos",
        default=2,
        type=check_positive_int,
        help="Threshold on the minimum number of read positions with at least"
        " 1 read for any junction to meet per-experiment filters for the LSVs"
        " it is a part of. Positions are relative to the aligned query"
        " sequences, ignoring soft clipping, and the first few bases of"
        " overhang on either end are ignored. When the minimum numbers of"
        " reads (--minreads) and positions are both met in enough experiments"
        " for a group in any junction that is part of an LSV, the LSV is"
        " considered admissible and associated coverage per experiment saved"
        " in the output MAJIQ files for potential downstream quantification."
        " [Default: %(default)s]",
    )

    # denovo flags
    buildparser_denovo = buildparser.add_argument_group("Denovo junctions options")
    buildparser_denovo.add_argument(
        "--min-denovo",
        default=5,
        type=check_positive_int,
        help="Threshold on the minimum total number of reads for a denovo"
        " junction to be detected for inclusion in the splicegraph. This"
        " per-experiment filter requires the --minpos filter to be satisfied"
        " at the same time. [Default: %(default)s]",
    )
    buildparser_denovo.add_argument(
        "--disable-denovo",
        dest="denovo",
        action="store_false",
        default=True,
        help="Disable denovo detection of junctions, splicesites and exons."
        " This will restrict analysis to junctions and exons found in provided"
        " annotations file, reducing the number of LSVs detected. Note that"
        " this does not disable detection of unannotated intron retention"
        " (see --disable-denovo-ir). [Default: denovo enabled]",
    )

    # intron retention flags
    buildparser_introns = buildparser.add_argument_group("Intron options")
    buildparser_introns.add_argument(
        "--irnbins",
        default=0.5,
        type=check_positive,
        help="Threshold on fraction of intronic read positions"
        " (aggregated/normalized to match junctions) with sufficient coverage"
        " (set by --min-intronic-cov) to pass per-experiment filters on"
        " introns. [Default: %(default)s]",
    )
    buildparser_introns.add_argument(
        "--min-intronic-cov",
        default=0.01,
        type=check_positive,
        help="Threshold on per-position normalized intronic readrate to be"
        " considered to have sufficient coverage at that position. Used with"
        " --irnbins to define per-experiment filters on introns."
        " [Default: %(default)s]",
    )
    buildparser_introns.add_argument(
        "--disable-ir",
        dest="ir",
        action="store_false",
        default=True,
        help="Disable intron retention detection. This applies to both"
        " annotated and unannotated retained introns."
        " [Default: intron retention enabled]",
    )
    buildparser_introns.add_argument(
        "--disable-denovo-ir",
        dest="denovo_ir",
        action="store_false",
        default=True,
        help="Disable detection of denovo introns only, keeping detection of"
        " annotated introns enabled. [Default: denovo introns enabled]",
    )
    buildparser_introns.add_argument(
        "--annotated_ir_always",
        dest="annot_ir_always",
        action="store_true",
        default=False,
        help="Automatically pass all annotated introns regardless of coverage."
        " By default, introns with insufficient coverage (not passing group"
        " filters) are excluded from the splicegraph and associated LSV"
        " definitions even if they are present in the annotations; this flag"
        " forces them to be kept.",
    )

    # incremental flags
    buildparser_incremental = buildparser.add_argument_group(
        "Incremental build options"
    )
    buildparser_incremental.add_argument(
        "--junc-files-only",
        dest="juncfiles_only",
        action="store_true",
        default=False,
        help="Stop MAJIQ builder execution after extracting junction"
        " information from BAM files into output SJ (*.sj) files for use in"
        " MAJIQ incremental builds. [Default: disabled]",
    )
    buildparser_incremental.add_argument(
        "--incremental",
        dest="aggregate",
        action="store_true",
        default=False,
        help="Enable use of SJ files generated by previous builds. This"
        " removes the need to reprocess the original BAM file again."
        " [Default: %(default)s]",
    )

    # simplifier flags
    buildparser_simplifier = buildparser.add_argument_group("Simplifier options")
    buildparser_simplifier.add_argument(
        "--simplify",
        dest="simpl_psi",
        default=-1,
        type=float,
        nargs="?",
        const=0.01,
        help="Enable the simplifier. Simplification ignores junctions and"
        " introns with consistently low usage within and between build groups"
        " in subsequent quantification. Optional threshold specifies maximum"
        " value of raw PSI considered as low usage (if not specified, uses"
        " %(const)s). If enabled, a junction/intron will be simplified if it"
        " has PSI above this threshold in less than min-experiments"
        " experiments for each build group and event it belongs to."
        " [Default: %(default)s]",
    )
    buildparser_simplifier.add_argument(
        "--simplify-min-experiments",
        default=None,
        type=check_positive,
        dest="simpl_min_exp",
        help="Override minimum number of experiments used for simplifier,"
        " alllowing for more relaxed or stringent number of experiments"
        " passing filters for/against simplification. Specified as with"
        " --min-experiments (fraction or absolute number). This is disabled by"
        " default and the value from --min-experiments is used instead.",
    )
    buildparser_simplifier.add_argument(
        "--simplify-annotated",
        dest="simpl_db",
        default=0,
        type=int,
        help="Simplifier minimum number or reads threshold for annotated"
        " junctions. If specified, an annotated junction must simultaneously"
        " have at least the specified number of reads (in addition to being"
        " above the PSI threshold from --simplify) in enough experiments to"
        " avoid removal by simplification. [Default: %(default)s]",
    )
    buildparser_simplifier.add_argument(
        "--simplify-denovo",
        dest="simpl_denovo",
        default=0,
        type=int,
        help="Simplifier minimum number or reads threshold for denovo"
        " junctions. If specified, a denovo junction must simultaneously have"
        "at least the specified number of reads (in addition to being above"
        " the PSI threshold from --simplify) in enough experiments to avoid"
        " removal by simplification. [Default: %(default)s]",
    )
    buildparser_simplifier.add_argument(
        "--simplify-ir",
        dest="simpl_ir",
        default=0,
        type=int,
        help="Simplifier minimum readrate threshold for intron retention."
        " If specified, a retained intron must simultaneously have at least"
        " the specified minimum normalized readrate (in addition to being"
        " above the PSI threshold from --simplify) in enough experiments to"
        " avoid removal by simplification. [Default: %(default)s]",
    )

    # bootstrapping
    buildparser_bootstrap = buildparser.add_argument_group(
        "Bootstrap coverage sampling"
    )
    buildparser_bootstrap.add_argument(
        "--markstacks",
        default=1e-7,
        type=float,
        dest="pvalue_limit",
        help="P-value threshold used for detecting and removing read stacks"
        " (outlier per-position read coverage under Poisson or"
        " negative-binomial null distribution). Use a negative value to"
        " disable stack detection/removal. [Default: %(default)s]",
    )
    buildparser_bootstrap.add_argument(
        "--m",
        default=30,
        type=check_positive_int,
        help="Number of bootstrap samples of total read coverage to save in"
        " output SJ and MAJIQ files for downstream quantification."
        " [Default: %(default)s]",
    )

    buildparser_advanced = buildparser.add_argument_group("Advanced options")

    # flag to save all possible LSVs
    buildparser_advanced.add_argument(
        "--permissive-lsvs",
        dest="lsv_strict",
        default=True,
        action="store_false",
        help="Consider all distinct LSVs, including events which are contained"
        " by other LSVs. By default, MAJIQ ignores all events for which their"
        " connections are all present in another event (def: redundant events)"
        " unless they are mutually redundant, in which case the events are"
        " equivalent (in this case the single-source event is selected)."
        " There are some cases where we would like to quantify these redundant"
        " events (excluding the equivalent mutually-redundant events); this"
        " flag enables more permissive output of splicing events.",
    )

    # flag to only output target LSVs
    buildparser_advanced.add_argument(
        "--target-lsvs",
        dest="only_target_lsvs",
        action="store_true",
        help="Only target LSVs will be to output majiq files",
    )

    # flag to only output source LSVs
    buildparser_advanced.add_argument(
        "--source-lsvs",
        dest="only_source_lsvs",
        action="store_true",
        help="Only source LSVs will be to output majiq files",
    )


    # flag to provide information about constitutive junctions
    buildparser_advanced.add_argument(
        "--dump-constitutive",
        dest="dump_const_j",
        action="store_true",
        default=False,
        help="Create constitutive_junctions.tsv file listing all junctions"
        " that pass group filters but are not part of any LSV because they"
        " are structurally constitutive. [Default: %(default)s]",
    )

    # flag to save per-position coverage for debugging to SJ files
    buildparser_advanced.add_argument(
        "--dump-coverage",
        dest="dump_coverage",
        action="store_true",
        default=False,
        help="Optionally dump raw junction coverage by position to created SJ"
        " files for experimental/debugging purposes",
    )

    sampling = new_subparser()

    sampling.add_argument(
        "--minreads",
        default=10,
        type=check_positive_int,
        help="Threshold on the minimum total number of reads for any junction"
        " or intron to meet per-experiment filters for the LSVs it is a part"
        " of. When the minimum numbers of reads and positions (--minpos) are"
        " both met in enough experiments for a group in any one"
        " junction/intron that is part of a LSV, the LSV is considered"
        " quantifiable in that group. [Default: %(default)s]",
    )
    sampling.add_argument(
        "--minpos",
        default=3,
        type=check_positive_int,
        help="Threshold on the minimum total number of read positions with at"
        " least 1 read for any junction or intron to meet per-experiment"
        " filters for the LSVs it is a part of. When the minimum number of"
        " reads (--minreads) and positions are both met in enough experiments"
        " for a group in any one junction/intron that is part of a LSV, the"
        " LSV is considered quantifiable in that group. [Default: %(default)s]",
    )

    psi = new_subparser()
    psi.add_argument(
        "files",
        nargs="+",
        action=StoreExperimentPaths,
        help="Paths to MAJIQ files for the experiment(s) to aggregate for"
        " PSI quantification as a single group.",
    )
    psi.add_argument(
        "-n",
        "--name",
        required=True,
        type=check_group_chars,
        help="Name used to identify the single group of experiments being"
        " quantified for output file name(s) and visualization using VOILA.",
    )

    psi.add_argument(
        "--output-type",
        choices=["voila", "tsv", "all"],
        default="all",
        help="Specify the type(s) of output files to produce: voila file to"
        " use with voila, TSV file with basic quantifications per LSV, or"
        " both. [Default: %(default)s]",
    )

    comparison = new_subparser()
    comparison_req = comparison.add_argument_group("Required specification of groups")
    comparison_req.add_argument(
        "-grp1",
        dest="files1",
        nargs="+",
        required=True,
        action=StoreExperimentPaths,
        help="Paths to MAJIQ files for the experiment(s) to quantify for first"
        " group (aggregated as replicates if deltapsi, independently if heterogen)",
    )
    comparison_req.add_argument(
        "-grp2",
        dest="files2",
        nargs="+",
        required=True,
        action=StoreExperimentPaths,
        help="Paths to MAJIQ files for the experiment(s) to quantify for first"
        " group (aggregated as replicates if deltapsi, independently if heterogen)",
    )
    comparison_req.add_argument(
        "-n",
        "--names",
        nargs=2,
        metavar=("NAME_GRP1", "NAME_GRP2"),
        required=True,
        action=StoreGroupNames,
        type=check_group_chars,
        help="The names that identify the groups being compared.",
    )

    delta = new_subparser()
    delta.add_argument(
        "--default-prior",
        action="store_true",
        default=False,
        help="Use a default prior instead of computing empirical prior from"
        " the data. [Default: default prior disabled]",
    )
    delta.add_argument(
        "--binsize",
        default=0.025,
        type=int,
        help="The bins for PSI values. With a --binsize of 0.025 (default), we"
        " have 40 bins. [Default: %(default)s]",
    )
    delta.add_argument(
        "--prior-minreads",
        default=20,
        type=check_positive_int,
        help="Minimum number of reads combining all positions in a junction to"
        " be considered (for the 'best set' calculation)."
        " [Default: %(default)s]",
    )
    delta.add_argument(
        "--prior-minnonzero",
        default=10,
        type=check_positive_int,
        help="Minimum number of positions for the best set. [Default: %(default)s]",
    )
    delta.add_argument(
        "--prior-iter",
        default=1,
        type=check_positive_int,
        dest="iter",
        help="Max number of iterations of the EM. [Default: %(default)s]",
    )
    delta.add_argument(
        "--output-type",
        choices=["voila", "tsv", "all"],
        default="all",
        help="Specify the type(s) of output files to produce: voila file to"
        " use with voila, TSV file with basic quantifications per LSV, or"
        " both. [Default: %(default)s]",
    )

    htrgen = new_subparser()
    htrgen.add_argument(
        "--keep-tmpfiles",
        action="store_true",
        default=False,
        dest="keep_tmpfiles",
        help="When this argument is specified, majiq heterogen will not remove"
        " the psi files that are temporary generated during the execution"
        " [Default: %(default)d]",
    )
    htrgen.add_argument(
        "--psi-samples",
        type=check_nonnegative_int,
        default=100,
        dest="psi_samples",
        help="Number of PSI samples to take per LSV junction. If equal to 0,"
        " use expected value only. [Default: %(default)d]",
    )
    htrgen.add_argument(
        "--stats",
        nargs="+",
        choices=["TTEST", "WILCOXON", "TNOM", "INFOSCORE", "ALL"],
        default=["TTEST", "WILCOXON", "TNOM"],
        help="Test statistics to run."
        " TTEST: unpaired two-sample t-test (Welch's t-test)."
        " WILCOXON: Mann-Whitney U two-sample test (nonparametric)."
        " TNOM: Total Number of Mistakes (nonparametric)."
        " INFOSCORE: TNOM but threshold maximizing mutual information with"
        " group labels (nonparametric)."
        " ALL: use all other available test statistics."
        " [Default: %(default)s]",
    )
    htrgen.add_argument(
        "--test_percentile",
        type=float,
        default=0.95,
        help="For each one of the statistical tests, we combine all pvalue per"
        " psi sample by percentile calculation. This argument allows the user"
        " define with percentile they want to use [Default: %(default)d]",
    )
    htrgen.add_argument(
        "--visualization-std",
        type=check_positive,
        default=1e-2,
        help="Change stochastic estimation error in terms of standard deviation"
        " of discretized average posterior per group by sampling additional"
        " values of PSI when number of samples is low [Default: %(default)f]",
    )

    # calcpsi flags
    subparsers = parser.add_subparsers(help="")
    parser_preprocess = subparsers.add_parser(
        "build",
        help="Preprocess SAM/BAM files as preparation for the rest of the"
        " tools (psi, deltapsi)",
        parents=[common, buildparser],
    )
    parser_preprocess.set_defaults(func=build)

    parser_calcpsi = subparsers.add_parser(
        "psi",
        help="Calculate aggregated PSI values for N experiments as a single"
        " group using MAJIQ files produced by `majiq build`",
        parents=[common, sampling, psi],
    )
    parser_calcpsi.set_defaults(func=calcpsi)

    parser_deltagroup = subparsers.add_parser(
        "deltapsi",
        help="Calculate Delta PSI values between two groups of experiments"
        " (treated as replicates) using MAJIQ files produced by `majiq build`",
        parents=[common, sampling, comparison, delta],
    )
    parser_deltagroup.set_defaults(func=deltapsi)

    parser_heterogen = subparsers.add_parser(
        "heterogen",
        help="Perform independent quantification and statistical comparison of"
        " PSI values between experiments/samples from two groups using MAJIQ"
        " files produced by `majiq build`",
        parents=[common, sampling, comparison, htrgen],
    )
    parser_heterogen.set_defaults(func=calc_independent)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    check_license(args.license)

    args.func(args)


if __name__ == "__main__":
    main()
