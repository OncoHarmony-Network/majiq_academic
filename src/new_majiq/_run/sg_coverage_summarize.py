"""
sg_coverage_summmarize.py

Summarize over sg_coverage

Author: Joseph K Aicher
"""

import argparse
import new_majiq as nm
import new_majiq.constants as constants

from new_majiq.experiments import bam_experiment_name
from new_majiq.logger import get_logger
from pathlib import Path
from new_majiq._run._majiq_args import check_nonnegative_factory
from new_majiq._run._run import GenericSubcommand
from typing import (
    List,
    Optional,
)


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
    if (missing := sorted(x for x in args.sg_coverage if not x.exists())) :
        raise ValueError(
            f"Was unable to find all input sg_coverage files ({missing = })"
        )
    if args.summary.exists():
        raise ValueError(f"Output {args.summary} already exists")

    log = get_logger()
    log.info("Lazily loading input sg_coverage")
    coverage = nm.MultiSpliceGraphReads(sorted(set(args.sg_coverage)))
    log.info(
        f"Taking the {args.reduction} over experiments for each"
        f" intron/junction and saving to {args.summary}"
    )
    summary = coverage.summarize_experiments(args.reduction, compute=False)
    (
        summary.expand_dims(experiment=[bam_experiment_name(args.summary)])
        .assign_attrs(
            bam_path=str(args.summary.resolve()),
            bam_version=nm.__version__,
            summarized_inputs=sorted(set(str(x.resolve()) for x in args.sg_coverage)),
            reduction=args.reduction,
        )
        .assign_coords(
            intron_hash=("experiment", [coverage.intron_checksum]),
            junction_hash=("experiment", [coverage.junction_checksum]),
        )
        .to_zarr(args.summary, mode="w", group=constants.NC_SGREADS)
    )

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
