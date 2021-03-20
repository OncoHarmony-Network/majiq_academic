"""
simplify.py

Update splicegraph with specified experiments split into simplifier groups

Author: Joseph K Aicher
"""

import argparse

import new_majiq.constants as constants

from new_majiq._run._majiq_args import check_nonnegative_factory
from pathlib import Path
from new_majiq._run._run import GenericSubcommand
from typing import (
    List,
    TYPE_CHECKING,
    Optional,
)

if TYPE_CHECKING:
    import pandas as pd


DESCRIPTION = "Update splicegraph with specified experiment groups"


def simplifier_threshold_args(
    parser: argparse.ArgumentParser,
    prefix: str = "",
) -> None:
    """arguments for simplifier thresholds

    Parameters
    ----------
    parser: argparse.ArgumentParser
        parser to add arguments to
    prefix: str
        add prefix to threshold command line arguments to avoid collisions
        (does not change destination variable)
    """
    thresholds = parser.add_argument_group("Simplifier filters")
    # min-experiments
    thresholds.add_argument(
        f"--{prefix}min-experiments",
        dest="simplify_min_experiments",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_SIMPLIFIER_MINEXPERIMENTS,
        help="Threshold for group filters. If < 1, the fraction of experiments"
        " in a group that must pass individual filters for a feature to be"
        " unsimplified. If greater, an absolute number. (default: %(default)s)",
    )
    # per-experiment thresholds
    thresholds.add_argument(
        f"--{prefix}minpsi",
        dest="simplify_minpsi",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_SIMPLIFIER_MINPSI,
        help="Minimum fraction of intron/junction readrates leaving or entering"
        " an exon in a single connection to count as evidence to unsimplify"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-annotated",
        dest="simplify_minreads_annotated",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_SIMPLIFIER_MINREADS_ANNOTATED,
        help="Minimum readrate for annotated junctions to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-denovo",
        dest="simplify_minreads_denovo",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_SIMPLIFIER_MINREADS_DENOVO,
        help="Minimum readrate for denovo junctions to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    thresholds.add_argument(
        f"--{prefix}minreads-ir",
        dest="simplify_minreads_ir",
        type=check_nonnegative_factory(float, False),
        default=constants.DEFAULT_SIMPLIFIER_MINREADS_INTRON,
        help="Minimum readrate for intron retention to count as evidence"
        " to unsimplify (default: %(default)s)",
    )
    return


def reset_simplified_args(parser: argparse.ArgumentParser) -> None:
    """do we reset simplified status in splicegraph?"""
    # do we reset simplified status?
    reset = parser.add_argument_group("Resetting connections to simplified")
    reset_ex = reset.add_mutually_exclusive_group()
    reset_ex.add_argument(
        "--reset-simplify",
        action="store_true",
        dest="reset_simplify",
        default=True,
        help="Reset all introns/junctions to simplified before processing input"
        " experiments (default: reset_simplify=%(default)s)",
    )
    reset_ex.add_argument(
        "--update-simplify",
        action="store_false",
        dest="reset_simplify",
        default=True,
        help="Continue unsimplifying from input splicegraph, do not reset"
        " simplification information (default: reset_simplify=%(default)s)",
    )
    return


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument("base_sg", type=Path, help="Path to base splicegraph")
    parser.add_argument(
        "grouped_experiments",
        type=Path,
        help="Path to TSV with required columns 'group' and 'sj' defining"
        " groups of experiments and the paths to their sj files. This does"
        " not have to be the same as what was used for building",
    )
    parser.add_argument("out_sg", type=Path, help="Path for output splicegraph")

    reset_simplified_args(parser)
    simplifier_threshold_args(parser)

    return


def get_grouped_experiments(path: Path) -> "pd.DataFrame":
    """Get grouped experiments from table

    Get grouped experiments from table. Verify that there are no duplicate
    experiments (insofar as their paths) and that the paths exist.
    """
    import pandas as pd

    df = pd.read_csv(path, sep="\t", usecols=["group", "sj"])
    # determine if any duplicated experiments
    duplicated_mask = df.sj.duplicated()
    if duplicated_mask.any():
        duplicated_sj = set(df.sj[duplicated_mask])
        raise ValueError(
            f"Requested simplification with repeated experiments {duplicated_sj}"
        )
    # verify that all paths exist
    for sj_path in df.sj:
        if not Path(sj_path).exists():
            raise ValueError(f"Unable to find input experiment {sj_path}")
    # return table of experiments/groups
    return df


def run(args: argparse.Namespace) -> None:
    if not args.base_sg.exists():
        raise ValueError(f"Unable to find base splicegraph at {args.base_sg}")
    if not args.grouped_experiments.exists():
        raise ValueError(
            f"Unable to find group definitions at {args.grouped_experiments}"
        )
    if args.out_sg.exists():
        raise ValueError(f"Output path {args.out_sg} already exists")
    # load experiments table (and verify correctness/paths exist)
    experiments = get_grouped_experiments(args.grouped_experiments)
    ngroups = experiments.group.nunique()

    # begin processing
    import new_majiq as nm
    from new_majiq.logger import get_logger

    log = get_logger()
    log.info(
        f"Updating {args.base_sg.resolve()} using"
        f" {len(experiments)} experiments in {ngroups} groups"
    )
    log.info("Loading base splicegraph")
    sg = nm.SpliceGraph.from_zarr(args.base_sg)

    if args.reset_simplify:
        log.info("Setting all introns and junctions to simplified")
        sg.introns._simplify_all()
        sg.junctions._simplify_all()

    log.info("Identifying introns and junctions to unsimplify")
    simplifier_group = sg.exon_connections.simplifier()
    for group_ndx, (group, group_sjs) in enumerate(experiments.groupby("group")["sj"]):
        log.info(
            f"Processing experiments from group {group} ({1 + group_ndx} / {ngroups})"
        )
        for sj_ndx, sj in enumerate(group_sjs):
            log.info(
                f"Processing coverage from {Path(sj).resolve()}"
                f" (experiment {1 + sj_ndx} / {len(group_sjs)} in group {group})"
            )
            simplifier_group.add_experiment(
                nm.SpliceGraphReads.from_connections_and_sj(
                    sg.introns,
                    sg.junctions,
                    nm.SJIntronsBins.from_zarr(sj),
                    nm.SJJunctionsBins.from_zarr(sj),
                ),
                min_psi=args.simplify_minpsi,
                minreads_annotated=args.simplify_minreads_annotated,
                minreads_denovo=args.simplify_minreads_denovo,
                minreads_introns=args.simplify_minreads_ir,
            )
        simplifier_group.update_connections(args.simplify_min_experiments)
    del simplifier_group

    log.info(f"Saving updated splicegraph to {args.out_sg.resolve()}")
    sg.to_zarr(args.out_sg)
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
