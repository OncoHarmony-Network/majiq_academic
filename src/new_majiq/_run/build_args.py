"""
build_args.py

Arguments for majiq-build components or pipeline

Author: Joseph K Aicher
"""

import argparse
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_nonnegative_factory,
)
from new_majiq.logger import get_logger

# type aliases
SJGroupT = Sequence[Path]
SJGroupsT = Dict[str, SJGroupT]


# what logic do we follow?
class BuildType(Enum):
    BUILD_ALL = "build_all"
    BUILD_KNOWN = "build_known"
    SIMPLIFY_ONLY = "simplify_only"


class IntronsType(Enum):
    NO_INTRONS = "no_introns"
    ANNOTATED_INTRONS = "annotated_introns"
    ALL_INTRONS = "all_introns"


def enable_simplify_args(parser: argparse.ArgumentParser) -> None:
    simplify = parser.add_argument_group("enable/disable simplifier arguments")
    simplify_ex = simplify.add_mutually_exclusive_group()
    simplify_ex.add_argument(
        "--simplify",
        dest="simplify",
        action="store_true",
        default=None,
        help="(Un)simplify splicegraph using evidence from input experiments"
        " (default: do not simplify unless build type set to --simplify-only)",
    )
    simplify_ex.add_argument(
        "--no-simplify",
        dest="simplify",
        action="store_false",
        default=None,
        help="Explicitly request to not perform simplification."
        " Will raise error if --simplify-only is set."
        " (default: do not simplify unless build type set to --simplify-only)",
    )
    return


def build_type_args(
    parser: argparse.ArgumentParser,
    allow_simplify_only: bool = True,
) -> None:
    """arguments for build type"""
    build = parser.add_argument_group("select build type arguments")
    build_ex = build.add_mutually_exclusive_group()
    build_ex.add_argument(
        "--build-all",
        dest="build",
        default=BuildType.BUILD_ALL,
        action="store_const",
        const=BuildType.BUILD_ALL,
        help="Process experiments to find new junctions and update existing junctions"
        " (default: %(default)s)",
    )
    build_ex.add_argument(
        "--build-known",
        dest="build",
        default=BuildType.BUILD_ALL,
        action="store_const",
        const=BuildType.BUILD_KNOWN,
        help="Process experiments to update known junctions"
        " (note that denovo/annotated intron processing specified separately)"
        " (default: %(default)s)",
    )
    if allow_simplify_only:
        build_ex.add_argument(
            "--simplify-only",
            dest="build",
            default=BuildType.BUILD_ALL,
            action="store_const",
            const=BuildType.SIMPLIFY_ONLY,
            help="Only perform simplification (default: %(default)s)",
        )
    return


def ir_filtering_args(parser: argparse.ArgumentParser) -> None:
    """arguments for intron filtering"""
    introns = parser.add_argument_group("intron filtering arguments")
    introns_ex = introns.add_mutually_exclusive_group()
    introns_ex.add_argument(
        "--all-introns",
        dest="introns",
        default=IntronsType.ALL_INTRONS,
        action="store_const",
        const=IntronsType.ALL_INTRONS,
        help="Keep all annotated introns and denovo introns passing build filters"
        " (default: %(default)s)",
    )
    introns_ex.add_argument(
        "--annotated-introns",
        dest="introns",
        default=IntronsType.ALL_INTRONS,
        action="store_const",
        const=IntronsType.ANNOTATED_INTRONS,
        help="Keep all annotated introns only (default: %(default)s)",
    )
    introns_ex.add_argument(
        "--no-introns",
        dest="introns",
        default=IntronsType.ALL_INTRONS,
        action="store_const",
        const=IntronsType.NO_INTRONS,
        help="Drop/ignore all introns (default: %(default)s)",
    )
    return


def build_threshold_args(parser: argparse.ArgumentParser) -> None:
    """arguments for build threshold parameters"""
    thresholds = parser.add_argument_group("build filters arguments")
    # min-experiments
    thresholds.add_argument(
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_BUILD_MINEXPERIMENTS,
        metavar="X",
        help="Threshold per experiments group. If < 1, the fraction of"
        " experiments in a group that must pass individual filters for a"
        " feature to be accepted. If greater, an absolute number."
        " (default: %(default)s)",
    )
    # per-experiment thresholds
    thresholds.add_argument(
        "--minreads",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_BUILD_MINREADS,
        metavar="N",
        help="Minimum readrate to pass an annotated junction (default: %(default)s)",
    )
    thresholds.add_argument(
        "--mindenovo",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_BUILD_MINDENOVO,
        metavar="N",
        help="Minimum readrate to pass a denovo junction or intron (default: %(default)s)",
    )
    thresholds.add_argument(
        "--minpos",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_BUILD_MINPOS,
        metavar="N",
        help="Minimum number of nonzero positions to pass a junction."
        " This is scaled for introns with some minimum coverage per bin to"
        " account for length dependence."
        " (default: %(default)s).",
    )
    thresholds.add_argument(
        "--max-pctbins",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_BUILD_MAX_PCTBINS,
        metavar="P",
        help="Maximum fraction of intron bins on which to require coverage"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        "--junction-acceptance-probability",
        metavar="P",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
        help="Set length-dependent minbins intron thresholds by considering"
        " per-position Poisson readrate for which junctions are accepted with"
        " this probability (default: %(default)s)",
    )
    thresholds.add_argument(
        "--intron-acceptance-probability",
        metavar="P",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_BUILD_MATCH_INTRON_PROBABILITY,
        help="Set length-dependent minbins intron thresholds by picking least"
        " strict thresholds for which per-position readrate determined by"
        " --junction-acceptance-probability gives has acceptance probability"
        " less than this probability (default: %(default)s)",
    )
    return


def reset_simplified_args(parser: argparse.ArgumentParser) -> None:
    """do we reset simplified status in splicegraph?"""
    # do we reset simplified status?
    reset = parser.add_argument_group("select if resetting simplification arguments")
    reset_ex = reset.add_mutually_exclusive_group()
    reset_ex.add_argument(
        "--reset-simplified",
        action="store_true",
        dest="reset_simplify",
        default=True,
        help="If simplifying, reset all introns/junctions to simplified,"
        " i.e. start simplification from scratch"
        " (default: reset_simplify=%(default)s)",
    )
    reset_ex.add_argument(
        "--update-simplified",
        action="store_false",
        dest="reset_simplify",
        default=True,
        help="Continue from existing simplified splicegraph,"
        " i.e. new denovo connections are simplified, existing connections will"
        " not be reset to simplified (default: reset_simplify=%(default)s)",
    )
    return


def simplifier_threshold_args(
    parser: argparse.ArgumentParser,
    prefix: str = "",
) -> None:
    """arguments for simplifier thresholds

    Parameters
    ----------
    parser: argparse.ArgumentParser
        parser to add arguments to
    prefix: str
        add prefix to threshold command line arguments to avoid collisions
        (does not change destination variable)
    """
    thresholds = parser.add_argument_group("simplifier filter arguments")
    # min-experiments
    thresholds.add_argument(
        f"--{prefix}min-experiments",
        metavar="X",
        dest="simplify_min_experiments",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_SIMPLIFIER_MINEXPERIMENTS,
        help="Threshold for group filters. If < 1, the fraction of experiments"
        " in a group that must pass individual filters for a feature to be"
        " unsimplified. If greater, an absolute number. (default: %(default)s)",
    )
    # per-experiment thresholds
    thresholds.add_argument(
        f"--{prefix}minpsi",
        metavar="P",
        dest="simplify_minpsi",
        type=check_nonnegative_factory(float, False),
        default=nm.constants.DEFAULT_SIMPLIFIER_MINPSI,
        help="Minimum fraction of intron/junction readrates leaving or entering"
        " an exon in a single connection to count as evidence to unsimplify"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-annotated",
        metavar="N",
        dest="simplify_minreads_annotated",
        type=check_nonnegative_factory(float, False),
        default=nm.constants.DEFAULT_SIMPLIFIER_MINREADS_ANNOTATED,
        help="Minimum readrate for annotated junctions to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-denovo",
        metavar="N",
        dest="simplify_minreads_denovo",
        type=check_nonnegative_factory(float, False),
        default=nm.constants.DEFAULT_SIMPLIFIER_MINREADS_DENOVO,
        help="Minimum readrate for denovo junctions to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-ir",
        metavar="N",
        dest="simplify_minreads_ir",
        type=check_nonnegative_factory(float, False),
        default=nm.constants.DEFAULT_SIMPLIFIER_MINREADS_INTRON,
        help="Minimum readrate for intron retention to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    return


def sj_strandness_args(parser: argparse.ArgumentParser) -> None:
    """arguments for configuring strandness for SJ parsing"""
    strandness = parser.add_argument_group("strandness arguments")
    strandness.add_argument(
        "--strandness",
        type=str,
        default="AUTO",
        choices=(
            "AUTO",
            nm.ExperimentStrandness.NONE.name,
            nm.ExperimentStrandness.FORWARD.name,
            nm.ExperimentStrandness.REVERSE.name,
        ),
        help="Strandness of input BAM."
        " AUTO = automatically detect strand (use median ratio of forward vs"
        " reverse stranded reads at annotated junctions). (default: %(default)s)",
    )
    strandness.add_argument(
        "--auto-minreads",
        metavar="N",
        type=check_nonnegative_factory(int, reject_zero=True),
        default=nm.constants.DEFAULT_BAM_STRAND_MINREADS,
        help="For automatic detection of strand. Only consider evidence from"
        " splicegraph junctions with at least this many total (unstranded) reads"
        " (default: %(default)s)",
    )
    strandness.add_argument(
        "--auto-minjunctions",
        metavar="N",
        type=check_nonnegative_factory(int, reject_zero=True),
        default=nm.constants.DEFAULT_BAM_STRAND_MINJUNCTIONS,
        help="For automatic detection of strand. Infer unstranded if the number"
        " of splicegraph junctions with sufficient reads is less than this argument"
        " (default: %(default)s)",
    )
    strandness.add_argument(
        "--auto-mediantolerance",
        metavar="X",
        type=check_nonnegative_factory(float, reject_zero=True),
        default=nm.constants.DEFAULT_BAM_STRAND_MINDEVIATION,
        help="For automatic detection of strand. Infer unstranded if the median"
        " proportion over junctions of forward strand reads vs all reads"
        " deviates from 0.5 by only this amount (default: %(default)s)",
    )
    return


def gff3_parsing_pipeline(
    gff3: Path,
    features_default: bool = True,
    features_tsv: Optional[Path] = None,
    types_genes: List[str] = [],
    types_transcripts: List[str] = [],
    types_exons: List[str] = [],
    types_silent: List[str] = [],
    types_hard_skip: List[str] = [],
) -> nm.SpliceGraph:
    """Parse GFF3 from specified path, return splicegraph"""
    log = get_logger()
    log.info("Loading annotated splicegraph from %s", gff3)

    def log_function(gff3_type: str, missing_reason: str, count: int) -> None:
        """log function for :meth:`SpliceGraph.from_gff3`

        log function for :meth:`SpliceGraph.from_gff3` as warnings for
        non-silent skipped parents/top-level ancestores
        """
        log.warning(
            "GFF3 type '%s' skipped %d times (%s) (see GFF3 parsing options)",
            gff3_type,
            count,
            missing_reason,
        )
        return

    # get GFF3 types for parsing
    gff3_types: nm.GFF3TypesMap
    if features_default and not features_tsv:
        gff3_types = nm.GFF3TypesMap()  # use default values
    else:
        gff3_types = nm.GFF3TypesMap({})  # start with empty values
        if features_tsv:
            # process each line in features_tsv
            for x in open(features_tsv, "r"):
                gff3_key, gff3_action, *_ = x.strip().split("\t")
                if gff3_key in gff3_types:
                    raise ValueError(
                        f"{features_tsv} specifies {gff3_key} more than once"
                    )
                gff3_types[gff3_key] = gff3_action
    gff3_types = nm.GFF3TypesMap.from_types_sets(
        gene_types=gff3_types.gene_types() | set(types_genes),
        transcript_types=gff3_types.transcript_types() | set(types_transcripts),
        exon_types=gff3_types.exon_types() | set(types_exons),
        silent_types=gff3_types.silent_types() | set(types_silent),
        hard_skip_types=gff3_types.hard_skip_types() | set(types_hard_skip),
    )
    log.debug("Parsing GFF3 with %s", gff3_types)

    sg = nm.SpliceGraph.from_gff3(
        gff3,
        process_ir=nm.constants.DEFAULT_BUILD_PROCESS_IR,
        gff3_types=gff3_types,
        log_function=log_function,
    )
    log.debug("Loaded %s", sg)
    return sg


def gff3_parsing_args(
    parser: argparse.ArgumentParser,
    prefix: str = "",
) -> None:
    """arguments for parsing GFF3 for transcript models

    Parameters
    ----------
    parser: argparse.ArgumentParser
        parser to add arguments to
    prefix: str
        add prefix to threshold command line arguments
        (does not change destination variable)
    """
    gff3_options = parser.add_argument_group("GFF3 parsing options arguments")
    gff3_base_types = gff3_options.add_mutually_exclusive_group()
    gff3_base_types.add_argument(
        f"--{prefix}features-default",
        dest="gff3_features_default",
        action="store_true",
        default=True,
        help="Initialize parsing of GFF3 transcripts using MAJIQ defaults"
        " for select GFF3 types."
        " Use --types-* to map additional GFF3 types. (default)",
    )
    gff3_base_types.add_argument(
        f"--{prefix}features-none",
        dest="gff3_features_default",
        action="store_false",
        default=True,
        help="Initialize parsing of GFF3 transcripts with no GFF3 types."
        " Use --types-* to map additional GFF3 types"
        " (default: features-default)",
    )
    gff3_base_types.add_argument(
        f"--{prefix}features-tsv",
        dest="gff3_features_tsv",
        metavar="tsv",
        default=None,
        type=ExistingResolvedPath,
        help="Initialize parsing of GFF3 transcripts with first two columns of"
        " tab delimited file (first column: GFF3 type to map, second column:"
        f" action for that GFF3 type, from {nm.GFF3TypesMap.VALID_ACTIONS()})"
        " (default: features-default)",
    )
    StoreGFF3Types = StoreRequiredUniqueActionFactory()
    gff3_options.add_argument(
        f"--{prefix}types-genes",
        dest="gff3_types_genes",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Map specified types from GFF3 as genes if found as ancestor of exon",
    )
    gff3_options.add_argument(
        f"--{prefix}types-transcripts",
        dest="gff3_types_transcripts",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Map specified types from GFF3 as transcripts if found as parent of exon",
    )
    gff3_options.add_argument(
        f"--{prefix}types-exons",
        dest="gff3_types_exons",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Map specified types from GFF3 as exons",
    )
    gff3_options.add_argument(
        f"--{prefix}types-silent",
        dest="gff3_types_silent",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Specified types from GFF3 will be ignored silently if found as"
        " parent (ignored transcript) or top-level ancestor (ignored gene)",
    )
    gff3_options.add_argument(
        f"--{prefix}types-hard-skip",
        dest="gff3_types_hard_skip",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Specified types from GFF3 should never be an ancestor of exons"
        " and can be completely skipped for potential performance improvements",
    )
    return


def _events_selection_args(
    parser: Union[argparse.ArgumentParser, argparse._ArgumentGroup]
) -> None:
    """selection of lsv events to produce output for"""
    select_lsvs = parser.add_mutually_exclusive_group()
    select_lsvs.add_argument(
        "--strict-lsvs",
        "--nonredundant-lsvs",
        dest="select_lsvs",
        default=nm.constants.DEFAULT_SELECT_LSVS,
        action="store_const",
        const=nm.constants.SelectLSVs.STRICT_LSVS,
        help="Select passed LSVs that are either not strict subsets of other"
        " events (nonredundant) or mutually redundant source events"
        " (i.e. strict LSVs) (default: %(default)s)",
    )
    select_lsvs.add_argument(
        "--permissive-lsvs",
        dest="select_lsvs",
        default=nm.constants.DEFAULT_SELECT_LSVS,
        action="store_const",
        const=nm.constants.SelectLSVs.PERMISSIVE_LSVS,
        help="Select all passed LSVs that are not mutually redundant targets"
        " (i.e. permissive LSVs) (default: %(default)s)",
    )
    select_lsvs.add_argument(
        "--source-lsvs",
        dest="select_lsvs",
        default=nm.constants.DEFAULT_SELECT_LSVS,
        action="store_const",
        const=nm.constants.SelectLSVs.SOURCE_LSVS,
        help="Select all passed LSVs that are source events (i.e. source LSVs)"
        " (default: %(default)s)",
    )
    select_lsvs.add_argument(
        "--target-lsvs",
        dest="select_lsvs",
        default=nm.constants.DEFAULT_SELECT_LSVS,
        action="store_const",
        const=nm.constants.SelectLSVs.TARGET_LSVS,
        help="Select all passed LSVs that are target events (i.e. target LSVs)"
        " (default: %(default)s)",
    )


def lsv_coverage_args(parser: argparse.ArgumentParser) -> None:
    """shared arguments for producing lsv coverage"""
    coverage = parser.add_argument_group("coverage arguments")
    coverage.add_argument(
        "--num-bootstraps",
        metavar="M",
        type=check_nonnegative_factory(int, False),
        default=nm.constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        help="Number of bootstrap replicates to sample (default: %(default)s)",
    )
    coverage.add_argument(
        "--stack-pvalue-threshold",
        metavar="P",
        type=check_nonnegative_factory(float, False),
        default=nm.constants.DEFAULT_COVERAGE_STACK_PVALUE,
        help="Exclude coverage from read stacks."
        " Read stacks are bins with readrate having right-tailed probability"
        " less than this threshold vs Poisson from other nonzero bins"
        " (default: %(default).2e)",
    )
    events = parser.add_argument_group("events selection arguments")
    events.add_argument(
        "--ignore-from",
        metavar="sg",
        type=ExistingResolvedPath,
        default=None,
        help="Path to other splicegraph, ignore LSVs shared with this splicegraph",
    )
    _events_selection_args(events)
    return


def do_build_junctions(
    sg: nm.SpliceGraph,
    experiments: SJGroupsT,
    experiment_thresholds: nm.ExperimentThresholds,
    min_experiments: float,
    process_denovo_junctions: bool,
    denovo_simplified: bool,
    imap_unordered_fn=map,
) -> nm.GeneJunctions:
    """Get updated GeneJunctions given input experiments and thresholds"""
    log = get_logger()
    log.info("Updating known %s and identifying denovo junctions", sg.junctions)
    junction_builder = sg.junctions.builder()
    for group_ndx, (group, group_sjs) in enumerate(experiments.items(), 1):
        log.info(
            "Processing junctions from group %s (%d / %d)",
            group,
            group_ndx,
            len(experiments),
        )
        build_group = sg.junctions.build_group(sg.exons)

        def add_experiment_to_group(sj: Path) -> Path:
            build_group.add_experiment(
                nm.SJJunctionsBins.from_zarr(sj),
                thresholds=experiment_thresholds,
                add_denovo=process_denovo_junctions,
            )
            return sj

        jobs = imap_unordered_fn(add_experiment_to_group, group_sjs)
        for sj_ndx, sj in enumerate(jobs, 1):
            log.info(
                "Processed junctions from %s (%d / %d)",
                sj,
                sj_ndx,
                len(group_sjs),
            )
        junction_builder.add_group(build_group, min_experiments)

    log.debug("Consolidating passed junctions from input experiments")
    result = junction_builder.get_passed(denovo_simplified)
    log.info("Updated %s from input experiments", result)
    return result


def do_build_introns(
    exons: nm.Exons,
    base_introns: nm.GeneIntrons,
    experiments: SJGroupsT,
    experiment_thresholds: nm.ExperimentThresholds,
    min_experiments: float,
    discard_denovo: bool,
    denovo_simplified: bool,
    imap_unordered_fn=map,
) -> nm.GeneIntrons:
    """Get updated GeneIntrons given input experiments and thresholds"""
    log = get_logger()
    log.info("Updating %s with %s and input coverage", base_introns, exons)
    log.debug("Updating flags with respect to annotated exon boundaries")
    base_introns = (
        exons.get_annotated()
        .potential_introns(make_simplified=denovo_simplified)
        .update_flags_from(base_introns)
    )
    log.debug("Identifying new passed introns")
    intron_group = base_introns.build_group()  # intron groups done in place

    def add_experiment_to_group(sj: Path) -> Path:
        intron_group.add_experiment(
            nm.SJIntronsBins.from_zarr(sj),
            thresholds=experiment_thresholds,
        )
        return sj

    for group_ndx, (group, group_sjs) in enumerate(experiments.items(), 1):
        log.info(
            "Processing introns from group %s (%d / %d)",
            group,
            group_ndx,
            len(experiments),
        )
        jobs = imap_unordered_fn(add_experiment_to_group, group_sjs)
        for sj_ndx, sj in enumerate(jobs, 1):
            log.info(
                "Processed introns from %s (%d / %d)",
                sj,
                sj_ndx,
                len(group_sjs),
            )
        intron_group.update_introns(min_experiments)
    log.debug("Propagating updated introns back to exons boundaries")
    result = base_introns.filter_passed(
        keep_annotated=True, discard_denovo=False
    ).propagate_to_annotated(
        annotated_exons=exons, keep_annotated=True, discard_denovo=discard_denovo
    )
    log.info("Updated %s from input experiments", result)
    return result


def do_build(
    sg: nm.SpliceGraph,
    experiments: SJGroupsT,
    experiment_thresholds: nm.ExperimentThresholds,
    min_experiments: float,
    process_denovo_junctions: bool,
    introns_type: IntronsType,
    denovo_simplified: bool,
    imap_unordered_fn=map,
) -> nm.SpliceGraph:
    """Get updated splicegraph given base splicegraph and experiments"""
    log = get_logger()

    updated_junctions = do_build_junctions(
        sg=sg,
        experiments=experiments,
        experiment_thresholds=experiment_thresholds,
        min_experiments=min_experiments,
        process_denovo_junctions=process_denovo_junctions,
        denovo_simplified=denovo_simplified,
        imap_unordered_fn=imap_unordered_fn,
    )

    log.info("Inferring denovo exons and updated exon boundaries")
    updated_exons = sg.exons.infer_with_junctions(updated_junctions)

    updated_introns: nm.GeneIntrons
    if introns_type == IntronsType.NO_INTRONS:
        log.debug("Dropping all introns")
        updated_introns = updated_exons.empty_introns()
    else:
        updated_introns = do_build_introns(
            exons=updated_exons,
            base_introns=sg.introns,
            experiments=experiments,
            experiment_thresholds=experiment_thresholds,
            min_experiments=min_experiments,
            discard_denovo=introns_type != IntronsType.ALL_INTRONS,
            denovo_simplified=denovo_simplified,
            imap_unordered_fn=imap_unordered_fn,
        )

    log.info("constructing splicegraph with updated exons and connections")
    return sg.with_updated_exon_connections(
        nm.ExonConnections.create_connecting(
            updated_exons, updated_introns, updated_junctions
        )
    )


def do_simplify(
    sg: nm.SpliceGraph,
    experiments: SJGroupsT,
    reset_simplify: bool,
    simplify_minpsi: float,
    simplify_minreads_annotated: float,
    simplify_minreads_denovo: float,
    simplify_minreads_ir: float,
    simplify_min_experiments: float,
    imap_unordered_fn=map,
) -> nm.SpliceGraph:
    """Update simplifier flags for input splicegraph in place given experiments"""
    log = get_logger()

    if reset_simplify:
        log.debug("Setting all introns and junctions to simplified")
        sg.introns._simplify_all()
        sg.junctions._simplify_all()

    simplifier_group = sg.exon_connections.simplifier()

    def add_experiment_to_group(sj: Path) -> Path:
        simplifier_group.add_experiment(
            nm.SJExperiment.from_zarr(sj),
            min_psi=simplify_minpsi,
            minreads_annotated=simplify_minreads_annotated,
            minreads_denovo=simplify_minreads_denovo,
            minreads_introns=simplify_minreads_ir,
        )
        return sj

    for group_ndx, (group, group_sjs) in enumerate(experiments.items(), 1):
        log.info(
            "Simplifying splicegraph with coverage from group %s (%d / %d)",
            group,
            group_ndx,
            len(experiments),
        )
        jobs = imap_unordered_fn(add_experiment_to_group, group_sjs)
        for sj_ndx, sj in enumerate(jobs, 1):
            log.info(
                "Simplified coverage with %s (%d / %d)",
                sj,
                sj_ndx,
                len(group_sjs),
            )
        simplifier_group.update_connections(simplify_min_experiments)

    return sg


def do_build_and_simplify(
    sg: nm.SpliceGraph,
    experiments: SJGroupsT,
    build: BuildType = BuildType.BUILD_ALL,
    minreads: int = nm.constants.DEFAULT_BUILD_MINREADS,
    mindenovo: int = nm.constants.DEFAULT_BUILD_DENOVO_JUNCTIONS,
    minpos: int = nm.constants.DEFAULT_BUILD_MINPOS,
    max_pctbins: float = nm.constants.DEFAULT_BUILD_MAX_PCTBINS,
    junction_acceptance_probability: float = nm.constants.DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
    intron_acceptance_probability: float = nm.constants.DEFAULT_BUILD_MATCH_INTRON_PROBABILITY,
    min_experiments: float = nm.constants.DEFAULT_BUILD_MINEXPERIMENTS,
    introns_type: IntronsType = IntronsType.ALL_INTRONS,
    simplify: bool = False,
    reset_simplify: bool = True,
    simplify_minpsi: float = nm.constants.DEFAULT_SIMPLIFIER_MINPSI,
    simplify_minreads_annotated: float = nm.constants.DEFAULT_SIMPLIFIER_MINREADS_ANNOTATED,
    simplify_minreads_denovo: float = nm.constants.DEFAULT_SIMPLIFIER_MINREADS_DENOVO,
    simplify_minreads_ir: float = nm.constants.DEFAULT_SIMPLIFIER_MINREADS_INTRON,
    simplify_min_experiments: float = nm.constants.DEFAULT_SIMPLIFIER_MINEXPERIMENTS,
    imap_unordered_fn=map,
) -> nm.SpliceGraph:
    log = get_logger()
    experiment_thresholds = nm.ExperimentThresholds(
        minreads=minreads,
        mindenovo=mindenovo,
        minpos=minpos,
        max_pctbins=max_pctbins,
        junction_acceptance_probability=junction_acceptance_probability,
        intron_acceptance_probability=intron_acceptance_probability,
    )
    if build:
        sg = do_build(
            sg=sg,
            experiments=experiments,
            experiment_thresholds=experiment_thresholds,
            min_experiments=min_experiments,
            # we are procesing denovo junction only if build_all specified
            process_denovo_junctions=build == BuildType.BUILD_ALL,
            introns_type=introns_type,
            # if simplifying or if --update-simplified (i.e. not --reset-simplified)
            denovo_simplified=simplify or not reset_simplify,
            imap_unordered_fn=imap_unordered_fn,
        )
    else:
        log.info("Filtering introns (%s)", introns_type)
        sg = nm.SpliceGraph.from_components(
            contigs=sg.contigs,
            genes=sg.genes,
            exons=sg.exons,
            junctions=sg.junctions,
            introns=(
                sg.introns.filter_passed(keep_annotated=True, discard_denovo=False)
                if introns_type == IntronsType.ALL_INTRONS
                else sg.introns.filter_passed(keep_annotated=True, discard_denovo=True)
                if introns_type == IntronsType.ANNOTATED_INTRONS
                else sg.exons.empty_introns()
            ),
        )

    # perform simplification?
    if simplify:
        sg = do_simplify(
            sg=sg,
            experiments=experiments,
            reset_simplify=reset_simplify,
            simplify_minpsi=simplify_minpsi,
            simplify_minreads_annotated=simplify_minreads_annotated,
            simplify_minreads_denovo=simplify_minreads_denovo,
            simplify_minreads_ir=simplify_minreads_ir,
            simplify_min_experiments=simplify_min_experiments,
            imap_unordered_fn=imap_unordered_fn,
        )

    # propagate intron flags within annotated exon boundaries
    if introns_type != IntronsType.NO_INTRONS:
        # propagate intron flags for introns within annotated exon boundaries
        sg = nm.SpliceGraph.from_components(
            contigs=sg.contigs,
            genes=sg.genes,
            exons=sg.exons,
            junctions=sg.junctions,
            introns=sg.introns.propagate_through_annotated(
                keep_annotated=True,
                discard_denovo=introns_type == IntronsType.ANNOTATED_INTRONS,
            ),
        )

    return sg


def quantifiability_threshold_args(
    parser: argparse.ArgumentParser,
    prefix: str = "",
) -> None:
    """arguments for quantifiability thresholds

    Parameters
    ----------
    parser: argparse.ArgumentParser
        parser to add arguments to
    prefix: str
        add prefix to threshold command line arguments to avoid collisions
        (does not change destination variable)
    """
    thresholds = parser.add_argument_group("quantifiability thresholds arguments")
    thresholds.add_argument(
        f"--{prefix}minreads",
        dest="quantify_minreads",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_QUANTIFY_MINREADS,
        help="Minimum readrate per experiment to pass a connection for quantification"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minbins",
        dest="quantify_minbins",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_QUANTIFY_MINBINS,
        help="Minimum number of nonzero bins to pass a connection for quantification"
        " (default: %(default)s).",
    )
