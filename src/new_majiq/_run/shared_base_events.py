"""
shared_base_events.py

Create zarr dataset with mask over e_idx and ec_idx for events found in base
splicegraph that are also found in derived splicegraphs

Author: Joseph K Aicher
"""

import argparse
import multiprocessing.dummy as mp
from typing import List, Optional

import dask.array as da
import numpy as np
import xarray as xr

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build_args import _events_selection_args
from new_majiq.logger import get_logger

DESCRIPTION = (
    "Create mask over e_idx and ec_idx for base events that are also found in"
    " derived splicegraphs"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "base_sg",
        type=ExistingResolvedPath,
        help="Path for base splicegraph introns, junctions",
    )
    parser.add_argument(
        "base_events",
        type=ExistingResolvedPath,
        help="Path with base events (PsiCoverage, PsiControlsSummary, etc.)",
    )
    parser.add_argument(
        "base_mask",
        type=NewResolvedPath,
        help="Path for mask over e_idx, ec_idx base events that are shared with"
        " derived splicegraphs",
    )
    parser.add_argument(
        "--sg-derived",
        metavar="SG",
        type=ExistingResolvedPath,
        required=True,
        nargs="+",
        help="Path for derived splicegraphs",
    )
    parser.add_argument(
        "--names-derived",
        metavar="name",
        type=str,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        required=True,
        help="Names matching derived splicegraph paths",
    )
    _events_selection_args(parser)
    resources_args(parser, use_dask=False)
    return


def run(args: argparse.Namespace) -> None:
    if len(args.sg_derived) != len(args.names_derived):
        raise ValueError(
            "--derived-sg and --derived-names must be given same number of arguments"
            f" ({len(args.sg_derived)} for --derived-sg,"
            f" {len(args.names_derived)} for --derived-names;"
            f" {args.sg_derived=}, {args.names_derived=})"
        )

    log = get_logger()

    log.info(
        f"Loading base events from {args.base_events} with splicegraph {args.base_sg}"
    )
    sg = nm.SpliceGraph.from_zarr(args.base_sg)
    events = nm.Events.from_zarr(args.base_events, sg.introns, sg.junctions)

    log.info(f"Initializing output mask over base events at {args.base_mask}")
    xr.Dataset(
        dict(
            ec_mask=(
                ("derived_name", "ec_idx"),
                da.empty(
                    (len(args.names_derived), events.num_connections),
                    dtype=bool,
                    chunks=(1, events.num_connections),
                ),
            ),
            e_mask=(
                ("derived_name", "e_idx"),
                da.empty(
                    (len(args.names_derived), events.num_events),
                    dtype=bool,
                    chunks=(1, events.num_events),
                ),
            ),
        ),
        dict(
            derived_name=args.names_derived,
            derived_paths=("derived_name", [str(x) for x in args.sg_derived]),
        ),
        dict(
            intron_hash=sg.introns.checksum(),
            junction_hash=sg.junctions.checksum(),
            derived_select=args.select_lsvs.name,
        ),
    ).to_zarr(args.base_mask, mode="w", consolidated=True, compute=False)

    def update_base_mask(idx: int):
        e_mask = np.zeros(events.num_events, dtype=bool)
        e_mask[
            nm.SpliceGraph.from_zarr(args.sg_derived[idx], genes=sg.genes)
            .exon_connections.lsvs(args.select_lsvs)
            .unique_events_mask(events)
            .shared_events_idx
        ] = True
        ec_mask = events.broadcast_eidx_to_ecidx(e_mask)
        return (
            xr.Dataset(
                dict(e_mask=("e_idx", e_mask), ec_mask=("ec_idx", ec_mask)),
            )
            .expand_dims(derived_name=1)
            .to_zarr(args.base_mask, region=dict(derived_name=slice(idx, 1 + idx)))
        )

    log.info(f"Saving masks for derived splicegraphs to {args.base_mask}")
    with mp.Pool(args.nthreads) as p:
        jobs = p.imap_unordered(update_base_mask, range(len(args.sg_derived)))
        for ndx, _ in enumerate(jobs, 1):
            log.info(
                f"Finished processing {ndx} / {len(args.sg_derived)}"
                " derived splicegraphs"
            )

    return


subcommand = GenericSubcommand(DESCRIPTION, add_args, run)


def main(sys_args: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    subcommand.add_args(parser)
    args = parser.parse_args(sys_args)
    subcommand.run(args)
    return


if __name__ == "__main__":
    main()
