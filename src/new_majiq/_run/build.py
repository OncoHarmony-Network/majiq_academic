"""
build.py

Update splicegraph with specified experiments split into build groups

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


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument("base_sg", type=Path, help="Path to base splicegraph")
    parser.add_argument(
        "grouped_experiments",
        type=Path,
        help="Path to TSV with required columns 'group' and 'sj' defining"
        " groups of experiments and the paths to their sj files",
    )
    parser.add_argument("out_sg", type=Path, help="Path for output splicegraph")

    # per-experiment thresholds
    parser.add_argument(
        "--minreads",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINREADS,
        help="Minimum readrate to pass an annotated junction (default: %(default)s)",
    )
    parser.add_argument(
        "--mindenovo",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINDENOVO,
        help="Minimum readrate to pass a denovo junction or intron (default: %(default)s)",
    )
    parser.add_argument(
        "--minpos",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINPOS,
        help="Minimum number of nonzero positions to pass a junction."
        " This is scaled for introns with some minimum coverage per bin to"
        " account for length dependence."
        " (default: %(default)s).",
    )
    parser.add_argument(
        "--max-pctbins",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MAX_PCTBINS,
        help="Maximum fraction of intron bins on which to require coverage"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--junction-acceptance-probability",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
        help="Set length-dependent minbins intron thresholds by considering"
        " per-position Poisson readrate for which junctions are accepted with"
        " this probability (default: %(default)s)",
    )
    parser.add_argument(
        "--intron-acceptance-probability",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MATCH_INTRON_PROBABILITY,
        help="Set length-dependent minbins intron thresholds by picking least"
        " strict thresholds for which per-position readrate determined by"
        " --junction-acceptance-probability gives has acceptance probability"
        " less than this probability (default: %(default)s)",
    )

    # min-experiments
    parser.add_argument(
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MINEXPERIMENTS,
        help="Threshold for group filters. If < 1, the fraction of experiments"
        " in a group that must pass individual filters for a feature to be"
        " accepted. If greater, an absolute number. (default: %(default)s)",
    )

    # denovo junctions
    denovo_junctions_ex = parser.add_mutually_exclusive_group()
    denovo_junctions_ex.add_argument(
        "--known-junctions-only",
        action="store_false",
        dest="process_denovo_junctions",
        default=constants.DEFAULT_BUILD_DENOVO_JUNCTIONS,
        help="Only process junctions already in base splicegraph"
        " (default: process_denovo_junctions=%(default)s)",
    )
    denovo_junctions_ex.add_argument(
        "--process-denovo-junctions",
        action="store_true",
        dest="process_denovo_junctions",
        default=constants.DEFAULT_BUILD_DENOVO_JUNCTIONS,
        help="Process all junctions, known and denovo."
        " (default: process_denovo_junctions=%(default)s)",
    )

    # denovo introns/annotated introns
    annotated_ir_ex = parser.add_mutually_exclusive_group()
    annotated_ir_ex.add_argument(
        "--keep-annotated-ir",
        action="store_true",
        dest="keep_annotated_ir",
        default=constants.DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        help="Keep annotated introns even if they did not have read support"
        " (default keep_annotated_ir = %(default)s)",
    )
    annotated_ir_ex.add_argument(
        "--only-passed-ir",
        action="store_false",
        dest="keep_annotated_ir",
        default=constants.DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        help="Only keep annotated introns if they pass build filters"
        " (default keep_annotated_ir = %(default)s)",
    )
    denovo_ir_ex = parser.add_mutually_exclusive_group()
    denovo_ir_ex.add_argument(
        "--process-denovo-introns",
        action="store_true",
        dest="process_denovo_introns",
        default=constants.DEFAULT_BUILD_DENOVO_IR,
        help="Process denovo introns that pass build filters"
        " (default process_denovo_introns = %(default)s)",
    )
    denovo_ir_ex.add_argument(
        "--ignore-denovo-introns",
        action="store_false",
        dest="process_denovo_introns",
        default=constants.DEFAULT_BUILD_DENOVO_IR,
        help="Ignore denovo introns regardless of evidence"
        " (default process_denovo_introns = %(default)s)",
    )

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
        raise ValueError(f"Requested build with repeated experiments {duplicated_sj}")
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
    from new_majiq.ExonConnections import ExonConnections
    from new_majiq.logger import get_logger

    log = get_logger()
    experiment_thresholds = nm.ExperimentThresholds(
        minreads=args.minreads,
        mindenovo=args.mindenovo,
        minpos=args.minpos,
        max_pctbins=args.max_pctbins,
        junction_acceptance_probability=args.junction_acceptance_probability,
        intron_acceptance_probability=args.intron_acceptance_probability,
    )
    log.info(
        f"Updating {args.base_sg.resolve()} using"
        f" {len(experiments)} experiments in {ngroups} groups"
    )
    log.info("Loading base splicegraph")
    sg = nm.SpliceGraph.from_zarr(args.base_sg)
    log.info("Updating known and identifying denovo junctions")
    junction_builder = sg.junctions.builder()
    for group_ndx, (group, group_sjs) in enumerate(experiments.groupby("group")["sj"]):
        log.info(f"Processing junctions from group {group} ({1 + group_ndx} / {ngroups})")
        build_group = sg.junctions.build_group(sg.exons)
        for sj_ndx, sj in enumerate(group_sjs):
            log.info(
                f"Processing junctions from {Path(sj).resolve()}"
                f" (experiment {1 + sj_ndx} / {len(group_sjs)} in group {group})"
            )
            build_group.add_experiment(
                nm.SJJunctionsBins.from_zarr(sj),
                thresholds=experiment_thresholds,
                add_denovo=args.process_denovo_junctions,
            )
        junction_builder.add_group(build_group, args.min_experiments)
    # get updated junctions
    log.info("Finalizing junctions accumulated across build groups")
    updated_junctions = junction_builder.get_passed()
    del build_group, junction_builder  # explicitly release these resources

    log.info("Inferring denovo exons and updated exon boundaries")
    updated_exons = sg.exons.infer_with_junctions(updated_junctions)

    log.info("Determining potential gene introns using updated exons")
    potential_introns = sg.introns.potential_introns(updated_exons)
    log.info("Identifying new passed introns")
    intron_group = potential_introns.build_group()  # intron groups done in place
    for group_ndx, (group, group_sjs) in enumerate(experiments.groupby("group")["sj"]):
        log.info(f"Processing introns from group {group} ({1 + group_ndx} / {ngroups})")
        for sj_ndx, sj in enumerate(group_sjs):
            log.info(
                f"Processing introns from {Path(sj).resolve()}"
                f" (experiment {1 + sj_ndx} / {len(group_sjs)} in group {group})"
            )
            intron_group.add_experiment(
                nm.SJIntronsBins.from_zarr(sj),
                thresholds=experiment_thresholds,
            )
        intron_group.update_introns(args.min_experiments)  # resets group too
    log.info("Filtering potential introns to those passing thresholds")
    updated_introns = potential_introns.filter_passed(
        keep_annotated=args.keep_annotated_ir,
        discard_denovo=not args.process_denovo_introns,
    )
    del intron_group, potential_introns  # explicitly release these resources

    log.info("Updating splicegraph with updated junctions, exons, and introns")
    sg = sg.with_updated_exon_connections(
        ExonConnections.create_connecting(updated_exons, updated_introns, updated_junctions)
    )

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
