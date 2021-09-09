"""
quantify.py

Quantify input event coverage files

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
        help="Paths to psi coverage files. All must have been generated"
        " using the same splicegraph",
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
        nargs="+",
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
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        help="Threshold for group filters. If < 1, the fraction of experiments"
        " in a group that must pass individual filters for a connection to be"
        " passed. If greater, an absolute number. (default: %(default)s)",
    )
    drop_ex = parser.add_mutually_exclusive_group()
    drop_ex.add_argument(
        "--drop-unquantifiable",
        action="store_true",
        dest="drop_unquantifiable",
        default=True,
        help="Drop records for unquantifiable events, saving storage space"
        " (default: drop_unquantifiable=%(default)s)",
    )
    drop_ex.add_argument(
        "--keep-unquantifiable",
        action="store_false",
        dest="drop_unquantifiable",
        default=True,
        help="Keep records for unquantifiable events, ensuring that files"
        " from different samples but sharing the same events are aligned"
        " (default: drop_unquantifiable=%(default)s)",
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

    log = get_logger()

    log.info(f"Opening coverage from {len(args.coverage)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.coverage)
    log.info(f"Aggregating coverage from {len(psicov.prefixes)} samples")
    psicov = psicov.sum("aggregate", min_experiments_f=args.min_experiments)
    if args.drop_unquantifiable:
        psicov = psicov.drop_unquantifiable()

    log.info("Quantifying coverage")
    q = psicov.quantifier_dataset(pmf_bins=args.pmf_bins, quantiles=args.quantiles)
    q.to_zarr(args.output, mode="w", group=constants.NC_EVENTSQUANTIFIED)
    psicov.events.to_zarr(args.output, mode="a", group=constants.NC_EVENTS)

    if args.tsv:
        log.info("Loading splicegraph for annotated TSV file")
        sg = nm.SpliceGraph.from_zarr(args.tsv[0])
        log.info("Saving annotated quantifications to TSV")
        psicov.as_dataframe(sg).to_csv(args.tsv[1], sep="\t")
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
