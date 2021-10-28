"""
sg_coverage_summmarize.py

Summarize over sg_coverage

Author: Joseph K Aicher
"""

import argparse
from pathlib import Path
from typing import List, Optional

import new_majiq as nm
import new_majiq.constants as constants
from new_majiq._run._majiq_args import check_nonnegative_factory
from new_majiq._run._run import GenericSubcommand
from new_majiq.experiments import bam_experiment_name
from new_majiq.logger import get_logger

DESCRIPTION = "Summarize splicegraph coverage from multiple experiments"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "summary",
        type=Path,
        help="Path for resulting summary splicegraph reads"
        " (summary experiment name will be inferred from the path)",
    )
    parser.add_argument(
        "sg_coverage",
        type=Path,
        nargs="+",
        help="Path for output coverage over introns/junctions",
    )
    parser.add_argument(
        "--reduction",
        type=str,
        choices=sorted({"sum", "mean", "median", "max"}),
        default="median",
        help="Summary function applied to coverage from input experiments (default: %(default)s)",
    )
    parser.add_argument(
        "--chunksize",
        type=check_nonnegative_factory(int, True),
        default=constants.NC_SGREADS_CHUNKS,
        help="Chunksize for summarized counts",
    )
    return


def run(args: argparse.Namespace) -> None:
    if missing := sorted(x for x in args.sg_coverage if not x.exists()):
        raise ValueError(
            f"Was unable to find all input sg_coverage files ({missing = })"
        )
    if args.summary.exists():
        raise ValueError(f"Output {args.summary} already exists")

    log = get_logger()
    log.info("Loading input splicegraph coverage")
    coverage = nm.SpliceGraphReads.from_zarr(sorted(set(args.sg_coverage)))
    log.info("Summarizing input splicegraph coverage")
    coverage = coverage.summarize(
        bam_experiment_name(args.summary), reduction=args.reduction
    )
    log.info(f"Saving summarized splicegraph coverage to {args.summary}")
    coverage.to_zarr(args.summary, mode="w", chunksize=args.chunksize)
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
