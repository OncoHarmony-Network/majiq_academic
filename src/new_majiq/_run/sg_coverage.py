"""
sg_coverage.py

Get coverage over splicegraph for a specific experiment

Author: Joseph K Aicher
"""

import argparse
from pathlib import Path
from typing import List, Optional

import new_majiq as nm
import new_majiq.constants as constants
from new_majiq._run._majiq_args import check_nonnegative_factory
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = (
    "Get raw coverage for input experiment at splicegraph introns and junctions,"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "sj",
        type=Path,
        help="Path to SJ coverage for experiment",
    )
    parser.add_argument(
        "splicegraph", type=Path, help="Path to splicegraph with introns/junctions"
    )
    parser.add_argument(
        "sg_coverage",
        type=Path,
        help="Path for output coverage over introns/junctions",
    )
    parser.add_argument(
        "--chunksize",
        type=check_nonnegative_factory(int, True),
        default=constants.NC_SGREADS_CHUNKS,
        help="Chunksize for per-experiment counts",
    )
    return


def run(args: argparse.Namespace) -> None:
    if not args.sj.exists():
        raise ValueError(f"Was unable to find input experiment at {args.sj}")
    if not args.splicegraph.exists():
        raise ValueError(f"Was unable to find splicegraph at {args.splicegraph}")
    if args.sg_coverage.exists():
        raise ValueError(f"Output path {args.sg_coverage} already exists")

    log = get_logger()
    log.info(f"Loading input splicegraph from {args.splicegraph.resolve()}")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    log.info(f"Loading input experiment from {args.sj.resolve()}")
    sj = nm.SJExperiment.from_zarr(args.sj)
    log.info("Obtaining coverage over introns and junctions")
    sg_coverage = nm.SpliceGraphReads.from_connections_and_sj(
        sg.introns, sg.junctions, sj
    )
    log.info(f"Saving coverage to {args.sg_coverage.resolve()}")
    sg_coverage.to_zarr(args.sg_coverage, mode="w", chunksize=args.chunksize)
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
