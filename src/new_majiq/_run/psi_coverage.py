"""
psi_coverage.py

Get coverage over LSVs for a specific experiment

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
    "Identify coverage for input experiment at splicegraph LSVs,"
    " removing per-bin stacks and bootstrapping replicates of coverage"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "splicegraph",
        type=Path,
        help="Path to splicegraph to define LSVs",
    )
    parser.add_argument(
        "psi_coverage",
        type=Path,
        help="Path for output psi-coverage file",
    )
    parser.add_argument(
        "sj",
        type=Path,
        nargs="+",
        help="Path to SJ coverage files for experiments",
    )
    thresholds = parser.add_argument_group("Thresholds for quantifiability")
    thresholds.add_argument(
        "--minreads",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINREADS,
        help="Minimum readrate per experiment to pass a connection"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        "--minbins",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINBINS,
        help="Minimum number of nonzero bins to pass a connection"
        " (default: %(default)s).",
    )
    coverage = parser.add_argument_group("Output coverage options")
    coverage.add_argument(
        "--num-bootstraps",
        type=check_nonnegative_factory(int, False),
        default=constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        help="Number of bootstrap replicates to sample (default: %(default)s)",
    )
    coverage.add_argument(
        "--stack-pvalue-threshold",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_COVERAGE_STACK_PVALUE,
        help="Bins with readrate having right-tailed probability less than this"
        " threshold vs Poisson from other nonzero bins will be ignored as"
        " outlier 'read stacks' (default: %(default).2e)",
    )
    events = parser.add_argument_group("Output event options")
    events.add_argument(
        "--ignore-from",
        metavar="sg",
        type=Path,
        default=None,
        help="Path to other splicegraph, ignore LSVs shared with this splicegraph",
    )
    return


def run(args: argparse.Namespace) -> None:
    if missing := sorted(set(x for x in args.sj if not x.exists())):
        raise ValueError(f"Unable to find all input SJ files ({missing =})")
    if len(unique := set(args.sj)) != len(args.sj):
        # get non-unique sj to report error
        non_unique = sorted(args.sj)
        for x in unique:
            non_unique.remove(x)  # removes first occurence
        raise ValueError(f"Non-unique input SJ files ({non_unique = })")
    if not args.splicegraph.exists():
        raise ValueError(f"Was unable to find splicegraph at {args.splicegraph}")
    if args.psi_coverage.exists():
        raise ValueError(f"Output path {args.psi_coverage} already exists")

    log = get_logger()
    log.info(f"Loading input splicegraph from {args.splicegraph.resolve()}")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    log.info("Defining LSVs for coverage")
    lsvs = sg.exon_connections.lsvs()
    if args.ignore_from is not None:
        log.info(f"Ignoring LSVs also found in {args.ignore_from.resolve()}")
        lsvs = lsvs[
            lsvs.unique_events_mask(
                nm.SpliceGraph.from_zarr(
                    args.ignore_from, genes=sg.genes
                ).exon_connections.lsvs()
            ).unique_events_mask
        ]
    for ndx, sj_path in enumerate(args.sj, 1):
        log.info(
            f"Loading SJ coverage from {sj_path.resolve()} ({ndx} / {len(args.sj)})"
        )
        sj_junctions = nm.SJJunctionsBins.from_zarr(sj_path)
        sj_introns = nm.SJIntronsBins.from_zarr(sj_path)
        log.info("Converting to LSV coverage")
        lsv_coverage = nm.EventsCoverage.from_events_and_sj(
            lsvs,
            sj_junctions,
            sj_introns,
            num_bootstraps=args.num_bootstraps,
            pvalue_threshold=args.stack_pvalue_threshold,
        )
        log.info("Converting to PSI coverage")
        psi_coverage = nm.PsiCoverage.from_events_coverage(
            lsv_coverage, args.minreads, args.minbins
        )
        log.info(f"Saving PSI coverage to {args.psi_coverage.resolve()}")
        psi_coverage.to_zarr(
            args.psi_coverage, append=(ndx > 1), consolidated=(ndx == len(args.sj))
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
