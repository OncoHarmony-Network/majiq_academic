"""
legacy_psi.py

Emulate the behavior of MAJIQ PSI from majiq v2.2

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

DESCRIPTION = (
    "Emulate legacy majiq psi: aggregated PSI values for a group of experiments"
    " to produce a voila file in the old format"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--nthreads", type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_QUANTIFY_NTHREADS,
        help="Number of threads to use (default: %(default)s)",
    )
    parser.add_argument(
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        help="Threshold for group filters. If < 1, the fraction of experiments"
        " in a group that must pass individual filters for a connection to be"
        " passed. If greater, an absolute number. (default: %(default)s)",
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
        "out_voila", type=Path,
        help="Path for output voila file",
    )
    parser.add_argument(
        "name", type=str,
        help="Name used to identify the single group of experiments being"
        " quantified for visualization with VOILA"
    )
    parser.add_argument(
        "splicegraph", type=Path,
        help="Path to splicegraph used to define events with coverage for"
        " annotating information required in the VOILA file",
    )
    parser.add_argument(
        "coverage", type=Path, nargs="+",
        help="Paths to events coverage files being quantified as a group",
    )
    return


def run(args: argparse.Namespace) -> None:
    import numpy as np
    import new_majiq as nm
    from new_majiq.logger import get_logger
    log = get_logger()

    try:
        if not all(p.exists() for p in args.coverage):
            missing = sorted(p for p in args.coverage if not p.exists())
            raise ValueError(f"Unable to find input coverage ({missing = })")
        if not args.splicegraph.exists():
            raise ValueError(f"Unable to find splicegraph file {args.splicegraph}")
        if args.out_voila.exists():
            raise ValueError(f"Output voila file {args.out_voila} already exists")
    except ValueError:
        log.exception("Error with specified input/output files:")
        raise

    try:
        from rna_voila.api import Matrix
        from rna_voila.constants import ANALYSIS_PSI, VOILA_FILE_VERSION
    except ModuleNotFoundError:
        log.exception("Saving VOILA file requires optional module rna_voila:")
        raise

    log.info(
        f"Quantifying group {args.name} to {args.out_voila} using splicegraph"
        f" {args.splicegraph} and {len(args.coverage)} coverage files"
    )
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
    log.info("Loading splicegraph to annotate quantified events")
    sg = nm.SpliceGraph.from_netcdf(args.splicegraph)
    q_events = q.get_events(sg.introns, sg.junctions)
    # get lsv id/description (legacy: lsv_type)
    event_id = sg.event_id(q_events.ref_exon_idx, q_events.event_type)
    event_description = sg.event_description(q_events.ref_exon_idx, q_events.event_type)
    # get start/end of each connection
    start = q_events.connection_start()
    end = q_events.connection_end()
    log.info("Performing quantifications")
    means = q.bootstrap_posterior_mean
    bins = q.bootstrap_discretized_pmf(nthreads=args.nthreads)
    # TODO get experiment names
    log.info(f"Saving quantifications to {args.out_voila}")
    with Matrix(args.out_voila, "w", voila_file=True, voila_tsv=False) as out_h5p:
        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_PSI
        out_h5p.experiment_names = [[b"TODO"]]
        out_h5p.group_names = [args.name]
        # each event
        for event_idx in range(q_events.num_events):
            connection_slice = q_events.connections_slice_for_event(event_idx)
            out_h5p.psi(event_id[event_idx]).add(
                lsv_type=event_description[event_idx],
                bins=bins[connection_slice],
                means=means[connection_slice],
                junctions=np.stack(
                    (start[connection_slice], end[connection_slice]),
                    axis=1,
                ),
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
