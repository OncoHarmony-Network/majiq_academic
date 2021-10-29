"""
psi_coverage.py

Get coverage over LSVs for a specific experiment

Author: Joseph K Aicher
"""

import argparse
import multiprocessing.dummy as mp
from pathlib import Path
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_nonnegative_factory,
    chunks_args,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.experiments import bam_experiment_name
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
    thresholds = parser.add_argument_group("quantifiability thresholds arguments")
    thresholds.add_argument(
        "--minreads",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_QUANTIFY_MINREADS,
        help="Minimum readrate per experiment to pass a connection"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        "--minbins",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_QUANTIFY_MINBINS,
        help="Minimum number of nonzero bins to pass a connection"
        " (default: %(default)s).",
    )
    coverage = parser.add_argument_group("coverage arguments")
    coverage.add_argument(
        "--num-bootstraps",
        metavar="M",
        type=check_nonnegative_factory(int, False),
        default=nm.constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        help="Number of bootstrap replicates to sample (default: %(default)s)",
    )
    coverage.add_argument(
        "--stack-pvalue-threshold",
        metavar="P",
        type=check_nonnegative_factory(float, False),
        default=nm.constants.DEFAULT_COVERAGE_STACK_PVALUE,
        help="Bins with readrate having right-tailed probability less than this"
        " threshold vs Poisson from other nonzero bins will be ignored as"
        " outlier 'read stacks' (default: %(default).2e)",
    )
    events = parser.add_argument_group("events selection arguments")
    events.add_argument(
        "--ignore-from",
        metavar="sg",
        type=ExistingResolvedPath,
        default=None,
        help="Path to other splicegraph, ignore LSVs shared with this splicegraph",
    )
    select_lsvs = events.add_mutually_exclusive_group()
    select_lsvs.add_argument(
        "--strict-lsvs",
        "--nonredundant-lsvs",
        dest="select_lsvs",
        default=nm.constants.DEFAULT_SELECT_LSVS,
        action="store_const",
        const=nm.constants.SelectLSVs.STRICT_LSVS,
        help="Select passed LSVs that are either not strict subsets of other"
        " events (nonredundant) or mutually redundant source events"
        " (i.e. strict LSVs) (default: %(default)s)",
    )
    select_lsvs.add_argument(
        "--permissive-lsvs",
        dest="select_lsvs",
        default=nm.constants.DEFAULT_SELECT_LSVS,
        action="store_const",
        const=nm.constants.SelectLSVs.PERMISSIVE_LSVS,
        help="Select all passed LSVs that are not mutually redundant targets"
        " (i.e. permissive LSVs) (default: %(default)s)",
    )
    select_lsvs.add_argument(
        "--source-lsvs",
        dest="select_lsvs",
        default=nm.constants.DEFAULT_SELECT_LSVS,
        action="store_const",
        const=nm.constants.SelectLSVs.SOURCE_LSVS,
        help="Select all passed LSVs that are source events (i.e. source LSVs)"
        " (default: %(default)s)",
    )
    select_lsvs.add_argument(
        "--target-lsvs",
        dest="select_lsvs",
        default=nm.constants.DEFAULT_SELECT_LSVS,
        action="store_const",
        const=nm.constants.SelectLSVs.TARGET_LSVS,
        help="Select all passed LSVs that are target events (i.e. target LSVs)"
        " (default: %(default)s)",
    )
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
                ).exon_connections.lsvs()
            ).unique_events_mask
        ]

    # define how to get from sj path to psicoverage file
    def sj_to_psicov(sj_path: Path) -> nm.PsiCoverage:
        return nm.PsiCoverage.from_sj_lsvs(
            nm.SJExperiment.from_zarr(sj_path),
            lsvs,
            minreads=args.minreads,
            minbins=args.minbins,
            num_bootstraps=args.num_bootstraps,
            pvalue_threshold=args.stack_pvalue_threshold,
        )

    if len(args.sj) == 1:
        # if there is only one file, don't bother with threads
        log.info(f"Inferring PSI coverage from {args.sj[0]}")
        psi_coverage = sj_to_psicov(args.sj[0])
        log.info(f"Saving PSI coverage to {args.psi_coverage}")
        psi_coverage.to_zarr(args.psi_coverage, ec_chunksize=args.chunksize)
    else:
        # precompute prefixes to use
        log.info("Precomputing prefixes corresponding to input SJ files")
        prefixes = [
            bam_experiment_name(nm.SJJunctionsBins.original_path_from_zarr(x))
            for x in args.sj
        ]
        # we have more than one input file
        log.info(f"Saving event information and metadata to {args.psi_coverage}")
        nm.PsiCoverage.to_zarr_slice_init(
            args.psi_coverage,
            lsvs.save_df,
            prefixes,
            args.num_bootstraps,
            psicov_attrs=dict(
                sj=[str(x) for x in args.sj],
                minreads=args.minreads,
                minbins=args.minbins,
            ),
            ec_chunksize=args.chunksize,
        )
        nm.rng_resize(args.nthreads)  # allow for as many rngs as there are threads
        with mp.Pool(args.nthreads) as p:
            jobs = p.imap_unordered(
                lambda x: (
                    sj_to_psicov(x[1]).to_zarr_slice(
                        args.psi_coverage,
                        slice(x[0], 1 + x[0]),
                        ec_chunksize=args.chunksize,
                    )
                ),
                list(enumerate(args.sj)),
            )
            for ndx, _ in enumerate(jobs, 1):
                log.info(f"Finished processing {ndx} / {len(args.sj)} SJ files")
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
