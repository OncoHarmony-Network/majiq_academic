"""
legacy_majiq.py

Output majiq file compatible with MAJIQ v2

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    lsv_coverage_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = "Output legacy lsv coverage file (.majiq) compatible with MAJIQ v2"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        help="Path to splicegraph to define LSVs",
    )
    parser.add_argument(
        "sj",
        type=ExistingResolvedPath,
        help="Path to SJ coverage file for output",
    )
    parser.add_argument(
        "majiq",
        type=NewResolvedPath,
        help="Path for new majiq file with LSV coverage",
    )
    lsv_coverage_args(parser)
    return


def run(args: argparse.Namespace) -> None:
    log = get_logger()
    log.info(f"Loading input splicegraph from {args.splicegraph}")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    log.info(f"Defining LSVs for coverage ({args.select_lsvs})")
    lsvs = sg.exon_connections.lsvs(args.select_lsvs)
    if args.ignore_from is not None:
        log.info(f"Ignoring LSVs also found in {args.ignore_from}")
        lsvs = lsvs[
            lsvs.unique_events_mask(
                nm.SpliceGraph.from_zarr(
                    args.ignore_from, genes=sg.genes
                ).exon_connections.lsvs(args.select_lsvs)
            ).unique_events_mask
        ]

    log.info("Loading and bootstrapping coverage over LSVs")
    coverage = nm.EventsCoverage.from_events_and_sj(
        lsvs,
        nm.SJExperiment.from_zarr(args.sj),
        num_bootstraps=args.num_bootstraps,
        pvalue_threshold=args.stack_pvalue_threshold,
    )
    log.info("Saving coverage to %s", args.majiq)
    coverage.to_legacy_majiq(args.majiq, exon_connections=sg.exon_connections)
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
