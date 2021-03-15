"""
lsv_coverage.py

Get coverage over LSVs for a specific experiment

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


DESCRIPTION = (
    "Identify coverage for input experiment at splicegraph LSVs,"
    " removing per-bin stacks and bootstrapping replicates of coverage"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "sj",
        type=Path,
        help="Path to SJ coverage for experiment",
    )
    parser.add_argument(
        "splicegraph",
        type=Path,
        help="Path to splicegraph to define LSVs",
    )
    parser.add_argument(
        "lsv_coverage",
        type=Path,
        help="Path for output LSV coverage files",
    )
    parser.add_argument(
        "--num-bootstraps",
        type=check_nonnegative_factory(int, False),
        default=constants.DEFAULT_BUILD_NUM_BOOTSTRAPS,
        help="Number of bootstrap replicates to sample (default: %(default)s)",
    )
    parser.add_argument(
        "--stack-pvalue-threshold",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_BUILD_STACK_PVALUE,
        help="Bins with readrate having right-tailed probability less than this"
        " threshold vs Poisson from other nonzero bins will be ignored as"
        " outlier 'read stacks' (default: %(default).2e)",
    )
    return


def run(args: argparse.Namespace) -> None:
    if not args.sj.exists():
        raise ValueError(f"Was unable to find input experiment at {args.sj}")
    if not args.splicegraph.exists():
        raise ValueError(f"Was unable to find splicegraph at {args.splicegraph}")
    if args.lsv_coverage.exists():
        raise ValueError(f"Output path {args.lsv_coverage} already exists")
    import new_majiq as nm
    from new_majiq.logger import get_logger

    log = get_logger()
    log.info(f"Loading input splicegraph from {args.splicegraph.resolve()}")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    log.info(f"Loading input experiment from {args.sj.resolve()}")
    sj_junctions = nm.SJJunctionsBins.from_zarr(args.sj)
    sj_introns = nm.SJIntronsBins.from_zarr(args.sj)
    log.info("Obtaining coverage over LSVs")
    lsv_coverage = nm.EventsCoverage.from_events_and_sj(
        sg.lsvs(),
        sj_junctions,
        sj_introns,
        num_bootstraps=args.num_bootstraps,
        pvalue_threshold=args.stack_pvalue_threshold,
    )
    log.info(f"Saving coverage to {args.lsv_coverage.resolve()}")
    lsv_coverage.to_zarr(args.lsv_coverage)
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
