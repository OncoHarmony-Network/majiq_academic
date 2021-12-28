"""
psi_coverage.py

Get coverage over LSVs for a specific experiment

Author: Joseph K Aicher
"""

import argparse
import multiprocessing.dummy as mp
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    chunks_args,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build_args import lsv_coverage_args, quantifiability_threshold_args
from new_majiq.logger import get_logger

DESCRIPTION = (
    "Identify coverage for input experiment at splicegraph LSVs,"
    " removing per-bin stacks and bootstrapping replicates of coverage"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        help="Path to splicegraph to define LSVs",
    )
    parser.add_argument(
        "psi_coverage",
        type=NewResolvedPath,
        help="Path for output psi-coverage file",
    )
    parser.add_argument(
        "sj",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Path to SJ coverage files for experiments",
    )
    quantifiability_threshold_args(parser)
    lsv_coverage_args(parser)
    chunks_args(parser, nm.constants.DEFAULT_COVERAGE_CHUNKS)
    resources_args(parser, use_dask=False)
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

    nm.rng_resize(args.nthreads)
    nm.PsiCoverage.convert_sj_batch(
        args.sj,
        lsvs,
        args.psi_coverage,
        minreads=args.quantify_minreads,
        minbins=args.quantify_minbins,
        num_bootstraps=args.num_bootstraps,
        pvalue_threshold=args.stack_pvalue_threshold,
        ec_chunksize=args.chunksize,
        imap_unordered_fn=map
        if len(args.sj) == 1
        else mp.Pool(args.nthreads).imap_unordered,
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
