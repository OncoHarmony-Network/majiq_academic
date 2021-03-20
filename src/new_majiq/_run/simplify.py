"""
simplify.py

Update splicegraph with specified experiments split into simplifier groups

Author: Joseph K Aicher
"""

import argparse

import new_majiq.constants as constants

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
        "--reset-simplify",
        action="store_true",
        dest="reset_simplify",
        default=True,
        help="Reset all introns/junctions to simplified before processing input"
        " experiments (default: reset_simplify=%(default)s)",
    )
    reset_ex.add_argument(
        "--update-simplify",
        action="store_false",
        dest="reset_simplify",
        default=True,
        help="Continue unsimplifying from input splicegraph, do not reset"
        " simplification information (default: reset_simplify=%(default)s)",
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

    reset_simplified_args(parser)
    simplifier_threshold_args(parser)

    return


def run(args: argparse.Namespace) -> None:
    if not args.base_sg.exists():
        raise ValueError(f"Unable to find base splicegraph at {args.base_sg}")
    if args.out_sg.exists():
        raise ValueError(f"Output path {args.out_sg} already exists")

    import new_majiq._run._build_pipeline as nm_build
    from new_majiq.logger import get_logger

    log = get_logger()
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

    sg = nm_build.simplify(
        sg=args.base_sg,
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
