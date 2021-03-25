"""
build.py

Update splicegraph with specified experiments split into build groups

Author: Joseph K Aicher
"""

import argparse

import new_majiq.constants as constants
import new_majiq as nm
import new_majiq._run._build_pipeline as nm_build

from new_majiq.logger import get_logger

from new_majiq._run._majiq_args import check_nonnegative_factory
from pathlib import Path
from new_majiq._run._run import GenericSubcommand
from typing import (
    List,
    Optional,
)


DESCRIPTION = "Update splicegraph with specified experiment groups"


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
    thresholds = parser.add_argument_group("Simplifier filters")
    # min-experiments
    thresholds.add_argument(
        f"--{prefix}min-experiments",
        dest="simplify_min_experiments",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_SIMPLIFIER_MINEXPERIMENTS,
        help="Threshold for group filters. If < 1, the fraction of experiments"
        " in a group that must pass individual filters for a feature to be"
        " unsimplified. If greater, an absolute number. (default: %(default)s)",
    )
    # per-experiment thresholds
    thresholds.add_argument(
        f"--{prefix}minpsi",
        dest="simplify_minpsi",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_SIMPLIFIER_MINPSI,
        help="Minimum fraction of intron/junction readrates leaving or entering"
        " an exon in a single connection to count as evidence to unsimplify"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-annotated",
        dest="simplify_minreads_annotated",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_SIMPLIFIER_MINREADS_ANNOTATED,
        help="Minimum readrate for annotated junctions to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-denovo",
        dest="simplify_minreads_denovo",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_SIMPLIFIER_MINREADS_DENOVO,
        help="Minimum readrate for denovo junctions to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-ir",
        dest="simplify_minreads_ir",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_SIMPLIFIER_MINREADS_INTRON,
        help="Minimum readrate for intron retention to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    return


def reset_simplified_args(parser: argparse.ArgumentParser) -> None:
    """do we reset simplified status in splicegraph?"""
    # do we reset simplified status?
    reset = parser.add_argument_group("Resetting connections to simplified")
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


def build_threshold_args(parser: argparse.ArgumentParser) -> None:
    """arguments for build threshold parameters"""
    thresholds = parser.add_argument_group("Build filters")
    # min-experiments
    thresholds.add_argument(
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MINEXPERIMENTS,
        help="Threshold per experiments group. If < 1, the fraction of"
        " experiments in a group that must pass individual filters for a"
        " feature to be accepted. If greater, an absolute number."
        " (default: %(default)s)",
    )
    # per-experiment thresholds
    thresholds.add_argument(
        "--minreads",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINREADS,
        help="Minimum readrate to pass an annotated junction (default: %(default)s)",
    )
    thresholds.add_argument(
        "--mindenovo",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINDENOVO,
        help="Minimum readrate to pass a denovo junction or intron (default: %(default)s)",
    )
    thresholds.add_argument(
        "--minpos",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINPOS,
        help="Minimum number of nonzero positions to pass a junction."
        " This is scaled for introns with some minimum coverage per bin to"
        " account for length dependence."
        " (default: %(default)s).",
    )
    thresholds.add_argument(
        "--max-pctbins",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MAX_PCTBINS,
        help="Maximum fraction of intron bins on which to require coverage"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        "--junction-acceptance-probability",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
        help="Set length-dependent minbins intron thresholds by considering"
        " per-position Poisson readrate for which junctions are accepted with"
        " this probability (default: %(default)s)",
    )
    thresholds.add_argument(
        "--intron-acceptance-probability",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MATCH_INTRON_PROBABILITY,
        help="Set length-dependent minbins intron thresholds by picking least"
        " strict thresholds for which per-position readrate determined by"
        " --junction-acceptance-probability gives has acceptance probability"
        " less than this probability (default: %(default)s)",
    )
    return


def ir_filtering_args(parser: argparse.ArgumentParser) -> None:
    """arguments for intron filtering"""
    introns_ex = parser.add_mutually_exclusive_group()
    introns_ex.add_argument(
        "--all-introns",
        dest="introns",
        default=nm_build.IntronsType.ALL_INTRONS,
        action="store_const",
        const=nm_build.IntronsType.ALL_INTRONS,
        help="Keep all annotated introns and denovo introns passing build filters"
        " (default: %(default)s)",
    )
    introns_ex.add_argument(
        "--annotated-introns",
        dest="introns",
        default=nm_build.IntronsType.ALL_INTRONS,
        action="store_const",
        const=nm_build.IntronsType.ANNOTATED_INTRONS,
        help="Keep all annotated introns only (default: %(default)s)",
    )
    introns_ex.add_argument(
        "--no-introns",
        dest="introns",
        default=nm_build.IntronsType.ALL_INTRONS,
        action="store_const",
        const=nm_build.IntronsType.NO_INTRONS,
        help="Drop/do not process all introns (default: %(default)s)",
    )
    return


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument("base_sg", type=Path, help="Path to base splicegraph")
    parser.add_argument("out_sg", type=Path, help="Path for output splicegraph")
    experiments_ex = parser.add_mutually_exclusive_group(required=True)
    experiments_ex.add_argument(
        "--groups-tsv",
        type=Path,
        help="Path to TSV with required columns 'group' and 'sj' defining"
        " groups of experiments and the paths to their sj files, allowing"
        " specification of multiple build groups simultaneously",
    )
    experiments_ex.add_argument(
        "--sjs",
        type=Path,
        nargs="+",
        help="Paths to input experiments as SJ files for a single build group",
    )

    build_ex = parser.add_mutually_exclusive_group()
    build_ex.add_argument(
        "--build-all",
        dest="build",
        default=nm_build.BuildType.BUILD_ALL,
        action="store_const",
        const=nm_build.BuildType.BUILD_ALL,
        help="Process experiments to find new junctions and update existing junctions"
        " (default: build=%(default)s)",
    )
    build_ex.add_argument(
        "--build-known",
        dest="build",
        default=nm_build.BuildType.BUILD_ALL,
        action="store_const",
        const=nm_build.BuildType.BUILD_KNOWN,
        help="Process experiments to update known junctions"
        " (note that denovo/annotated intron processing specified separately)"
        " (default: build=%(default)s)",
    )
    build_ex.add_argument(
        "--simplify-only",
        dest="build",
        default=nm_build.BuildType.BUILD_ALL,
        action="store_const",
        const=nm_build.BuildType.BUILD_KNOWN,
        help="Only perform simplification (default: build=%(default)s)",
    )
    ir_filtering_args(parser)
    parser.add_argument(
        "--simplify",
        action="store_true",
        default=False,
        help="(Un)simplify splicegraph using evidence from input experiments"
        " (default: %(default)s)",
    )

    reset_simplified_args(parser)
    build_threshold_args(parser)
    simplifier_threshold_args(parser, prefix="simplify-")
    return


def run(args: argparse.Namespace) -> None:
    if not args.base_sg.exists():
        raise ValueError(f"Unable to find base splicegraph at {args.base_sg}")
    if args.out_sg.exists():
        raise ValueError(f"Output path {args.out_sg} already exists")

    log = get_logger()
    experiment_thresholds = nm.ExperimentThresholds(
        minreads=args.minreads,
        mindenovo=args.mindenovo,
        minpos=args.minpos,
        max_pctbins=args.max_pctbins,
        junction_acceptance_probability=args.junction_acceptance_probability,
        intron_acceptance_probability=args.intron_acceptance_probability,
    )
    if args.groups_tsv:
        if not args.groups_tsv.exists():
            raise ValueError(
                f"Unable to find group definitions at {args.groups_tsv.resolve()}"
            )
        log.info("Loading experiment groups")
        experiments = nm_build.get_grouped_experiments(args.groups_tsv)
    else:
        if (missing := sorted(set(x for x in args.sjs if not x.exists()))) :
            raise ValueError(f"Unable to find all input SJ files ({missing =})")
        experiments = {"": args.sjs}
    log.info(f"Loading base splicegraph from {args.base_sg.resolve()}")
    sg = nm.SpliceGraph.from_zarr(args.base_sg)

    # determine if we are simplifying in the end
    simplify = args.simplify or args.build == nm_build.BuildType.SIMPLIFY_ONLY

    # perform build?
    if args.build != nm_build.BuildType.SIMPLIFY_ONLY:
        sg = nm_build.build(
            sg=sg,
            experiments=experiments,
            experiment_thresholds=experiment_thresholds,
            min_experiments=args.min_experiments,
            process_denovo_junctions=args.build == nm_build.BuildType.BUILD_ALL,
            introns_type=args.introns,
            # if simplifying or if --update-simplified (i.e. not --reset-simplified)
            denovo_simplified=simplify or not args.reset_simplify,
        )
    else:
        log.info(f"Filtering introns {args.introns = }")
        sg = nm.SpliceGraph.from_components(
            contigs=sg.contigs,
            genes=sg.genes,
            exons=sg.exons,
            junctions=sg.junctions,
            introns=(
                sg.introns.filter_passed(keep_annotated=True, discard_denovo=False)
                if args.introns == nm_build.IntronsType.ALL_INTRONS else
                sg.introns.filter_passed(keep_annotated=True, discard_denovo=True)
                if args.introns == nm_build.IntronsType.ANNOTATED_INTRONS else
                sg.exons.empty_introns()
            )
        )

    # perform simplification?
    if simplify:
        sg = nm_build.simplify(
            sg=sg,
            experiments=experiments,
            reset_simplify=args.reset_simplify,
            simplify_minpsi=args.simplify_minpsi,
            simplify_minreads_annotated=args.simplify_minreads_annotated,
            simplify_minreads_denovo=args.simplify_minreads_denovo,
            simplify_minreads_ir=args.simplify_minreads_ir,
            simplify_min_experiments=args.simplify_min_experiments,
        )

    log.info(f"Saving updated splicegraph to {args.out_sg.resolve()}")
    sg.to_zarr(args.out_sg)
    return


def main(sys_args: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args(sys_args)
    run(args)
    return


subcommand = GenericSubcommand(DESCRIPTION, add_args, run)


if __name__ == "__main__":
    main()
