"""
build.py

Update splicegraph with specified experiments split into build groups

Author: Joseph K Aicher
"""

import argparse

import new_majiq.constants as constants

from new_majiq._run._majiq_args import check_nonnegative_factory
from pathlib import Path
from new_majiq._run._run import GenericSubcommand
from typing import (
    List,
    TYPE_CHECKING,
    Optional,
)

if TYPE_CHECKING:
    import pandas as pd


DESCRIPTION = "Update splicegraph with specified experiment groups"


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
    """add arguments to parser for filtering annotated/denovo introns"""
    ir_filtering = parser.add_argument_group("Intron retention filtering")
    # denovo introns/annotated introns
    annotated_ir_ex = ir_filtering.add_mutually_exclusive_group()
    annotated_ir_ex.add_argument(
        "--keep-annotated-ir",
        action="store_true",
        dest="keep_annotated_ir",
        default=constants.DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        help="Keep annotated introns even if they did not have read support"
        " (default keep_annotated_ir = %(default)s)",
    )
    annotated_ir_ex.add_argument(
        "--only-passed-ir",
        action="store_false",
        dest="keep_annotated_ir",
        default=constants.DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        help="Only keep annotated introns if they pass build filters"
        " (default keep_annotated_ir = %(default)s)",
    )
    denovo_ir_ex = ir_filtering.add_mutually_exclusive_group()
    denovo_ir_ex.add_argument(
        "--process-denovo-introns",
        action="store_true",
        dest="process_denovo_introns",
        default=constants.DEFAULT_BUILD_DENOVO_IR,
        help="Process denovo introns that pass build filters"
        " (default process_denovo_introns = %(default)s)",
    )
    denovo_ir_ex.add_argument(
        "--ignore-denovo-introns",
        action="store_false",
        dest="process_denovo_introns",
        default=constants.DEFAULT_BUILD_DENOVO_IR,
        help="Ignore denovo introns regardless of evidence"
        " (default process_denovo_introns = %(default)s)",
    )
    return


def denovo_junctions_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for processing of denovo junctions"""
    denovo_junctions = parser.add_argument_group("Denovo junction filtering")
    # denovo junctions
    denovo_junctions_ex = denovo_junctions.add_mutually_exclusive_group()
    denovo_junctions_ex.add_argument(
        "--known-junctions-only",
        action="store_false",
        dest="process_denovo_junctions",
        default=constants.DEFAULT_BUILD_DENOVO_JUNCTIONS,
        help="Only process junctions already in base splicegraph"
        " (default: process_denovo_junctions=%(default)s)",
    )
    denovo_junctions_ex.add_argument(
        "--process-denovo-junctions",
        action="store_true",
        dest="process_denovo_junctions",
        default=constants.DEFAULT_BUILD_DENOVO_JUNCTIONS,
        help="Process all junctions, known and denovo."
        " (default: process_denovo_junctions=%(default)s)",
    )
    return


def denovo_simplified_args(parser: argparse.ArgumentParser) -> None:
    """arguments for if denovos added simplified or not"""
    denovo_simplified = parser.add_argument_group("Simplification of new denovos")
    denovo_simplified_ex = denovo_simplified.add_mutually_exclusive_group()
    denovo_simplified_ex.add_argument(
        "--denovo-simplified",
        action="store_true",
        dest="denovo_simplified",
        default=constants.DEFAULT_BUILD_DENOVO_SIMPLIFIED,
        help="Denovo introns/junctions will be initially added marked as"
        " simplified (default: denovo_simplified=%(default)s)",
    )
    denovo_simplified_ex.add_argument(
        "--denovo-unsimplified",
        action="store_true",
        dest="denovo_simplified",
        default=constants.DEFAULT_BUILD_DENOVO_SIMPLIFIED,
        help="Denovo introns/junctions will be initially added as unsimplified"
        " (default: denovo_simplified=%(default)s)",
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

    build_threshold_args(parser)
    denovo_junctions_args(parser)
    ir_filtering_args(parser)
    denovo_simplified_args(parser)
    return


def run(args: argparse.Namespace) -> None:
    if not args.base_sg.exists():
        raise ValueError(f"Unable to find base splicegraph at {args.base_sg}")
    if args.out_sg.exists():
        raise ValueError(f"Output path {args.out_sg} already exists")

    import new_majiq as nm
    import new_majiq._run._build_pipeline as nm_build
    from new_majiq.logger import get_logger

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
        if (missing := sorted(set(x for x in args.sjs if not x.exists()))):
            raise ValueError(f"Unable to find all input SJ files ({missing =})")
        experiments = {"": args.sjs}

    sg = nm_build.build(
        sg=args.base_sg,
        experiments=experiments,
        experiment_thresholds=experiment_thresholds,
        min_experiments=args.min_experiments,
        process_denovo_junctions=args.process_denovo_junctions,
        process_denovo_introns=args.process_denovo_introns,
        keep_annotated_ir=args.keep_annotated_ir,
        denovo_simplified=args.denovo_simplified,
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
