"""
quantify.py

Quantify input event coverage files

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

DESCRIPTION = "Quantify events coverage files"


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "output",
        type=Path,
        help="Path for new-majiq quantification format",
    )
    parser.add_argument(
        "coverage",
        type=Path,
        nargs="+",
        help="Paths to events coverage files. All must have been generated"
        " using the same splicegraph",
    )
    parser.add_argument(
        "--nthreads",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_QUANTIFY_NTHREADS,
        help="Number of threads used for pmf bins and quantiles,"
        " which require more computation (default: %(default)s)",
    )
    parser.add_argument(
        "--pmf-bins",
        type=check_nonnegative_factory(int, True),
        default=None,
        help="If specified, number of discretized probability bins to save"
        " posterior distribution to (default: %(default)s)",
    )
    parser.add_argument(
        "--quantiles",
        type=float,
        default=[],
        help="Posterior quantiles to calculate (default: %(default)s)",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=None,
        nargs=2,
        metavar=("splicegraph", "tsv"),
        help="If specified, save quantifications to TSV annotated by input"
        " splicegraph",
    )
    parser.add_argument(
        "--minreads",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINREADS,
        help="Minimum readrate per experiment to pass a connection"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--minbins",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINBINS,
        help="Minimum number of nonzero bins to pass a connection"
        " (default: %(default)s).",
    )
    parser.add_argument(
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        help="Threshold for group filters. If < 1, the fraction of experiments"
        " in a group that must pass individual filters for a connection to be"
        " passed. If greater, an absolute number. (default: %(default)s)",
    )
    return


def run(args: argparse.Namespace) -> None:
    if not all(p.exists() for p in args.coverage):
        missing = sorted(p for p in args.coverage if not p.exists())
        raise ValueError(f"Unable to find input coverage ({missing = })")
    if args.output.exists():
        raise ValueError(f"Output coverage file {args.output} already exists")
    if args.tsv:
        if not args.tsv[0].exists():  # splicegraph
            raise ValueError(f"Unable to find splicegraph {args.tsv[0]} for TSV output")
        if args.tsv[1].exists():
            raise ValueError(f"Output TSV file {args.tsv[1]} already exists")

    import new_majiq as nm
    from new_majiq.logger import get_logger

    log = get_logger()

    thresholds = nm.QuantifierThresholds(
        minreads=args.minreads,
        minbins=args.minbins,
        min_experiments_f=args.min_experiments,
    )
    log.info("Determining quantifiable events")
    quantifiable = nm.QuantifiableEvents.from_quantifier_group(
        args.coverage, thresholds=thresholds
    )
    log.info("Aggregating coverage at these events for quantification")
    q = nm.QuantifiableCoverage.from_quantifier_group(
        args.coverage, quantifiable=quantifiable
    )
    log.info("Quantifying using input coverage and saving to file")
    q.to_zarr(
        args.output,
        pmf_bins=args.pmf_bins,
        quantiles=args.quantiles,
        nthreads=args.nthreads,
    )
    if args.tsv:
        log.info("Loading splicegraph for annotated TSV file")
        sg = nm.SpliceGraph.from_zarr(args.tsv[0])
        log.info("Saving annotated quantifications to TSV")
        q.as_dataframe(sg).to_csv(args.tsv[1], sep="\t")
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
