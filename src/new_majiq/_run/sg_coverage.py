"""
sg_coverage.py

Get coverage over splicegraph for a specific experiment

Author: Joseph K Aicher
"""

import argparse
import multiprocessing.dummy as mp
from pathlib import Path
from typing import List, Optional

import new_majiq as nm
import new_majiq.constants as constants
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_nonnegative_factory,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.experiments import bam_experiment_name
from new_majiq.logger import get_logger

DESCRIPTION = (
    "Get raw coverage for input experiment at splicegraph introns and junctions,"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        help="Path to splicegraph with introns/junctions",
    )
    parser.add_argument(
        "sg_coverage",
        type=NewResolvedPath,
        help="Path for output coverage over introns/junctions",
    )
    parser.add_argument(
        "sj",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Path to SJ coverage for experiments",
    )
    parser.add_argument(
        "--chunksize",
        type=check_nonnegative_factory(int, True),
        default=constants.NC_SGREADS_CHUNKS,
        help="Chunksize for per-experiment counts",
    )
    parser.add_argument(
        "--nthreads",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BAM_NTHREADS,
        help="Number of threads used for simultaneous processing of multiple"
        " input SJ files (default: %(default)s)",
    )
    return


def run(args: argparse.Namespace) -> None:
    log = get_logger()
    log.info(f"Loading input splicegraph from {args.splicegraph}")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)

    def sj_to_sgreads(sj_path: Path) -> nm.SpliceGraphReads:
        return nm.SpliceGraphReads.from_connections_and_sj(
            sg.introns, sg.junctions, nm.SJExperiment.from_zarr(sj_path)
        )

    if len(args.sj) == 1:
        log.info(f"Inferring SpliceGraphReads from {args.sj[0]}")
        sgreads = sj_to_sgreads(args.sj[0])
        log.info(f"Saving SpliceGraphReads to {args.sg_coverage}")
        sgreads.to_zarr(args.sg_coverage, chunksize=args.chunksize)
    else:
        # precompute prefixes to use
        log.info("Precomputing prefixes corresponding to input SJ files")
        prefixes = [
            bam_experiment_name(nm.SJJunctionsBins.original_path_from_zarr(x))
            for x in args.sj
        ]
        log.info(f"Saving prefixes and metadata to {args.sg_coverage}")
        nm.SpliceGraphReads.to_zarr_slice_init(
            args.sg_coverage,
            prefixes,
            len(sg.introns),
            len(sg.junctions),
            chunksize=args.chunksize,
            attrs=dict(
                sj=[str(x) for x in args.sj],
                sg=str(args.splicegraph),
            ),
        )
        with mp.Pool(args.nthreads) as p:
            jobs = p.imap_unordered(
                lambda x: (
                    sj_to_sgreads(x[1]).to_zarr_slice(
                        args.sg_coverage,
                        slice(x[0], 1 + x[0]),
                        chunksize=args.chunksize,
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
