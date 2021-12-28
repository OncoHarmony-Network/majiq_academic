"""
sg_coverage.py

Get coverage over splicegraph for a specific experiment

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
from new_majiq.logger import get_logger

DESCRIPTION = (
    "Get raw coverage for input experiment at splicegraph introns and junctions"
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
    chunks_args(parser, nm.constants.NC_SGREADS_CHUNKS)
    resources_args(parser, use_dask=False)
    return


def run(args: argparse.Namespace) -> None:
    log = get_logger()
    log.info(f"Loading input splicegraph from {args.splicegraph}")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    nm.SpliceGraphReads.convert_sj_batch(
        args.sj,
        sg.introns,
        sg.junctions,
        args.sg_coverage,
        chunksize=args.chunksize,
        attrs=dict(sg=str(args.splicegraph)),
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
