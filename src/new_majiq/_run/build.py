"""
build.py

Update splicegraph with specified experiments split into build groups

Author: Joseph K Aicher
"""

import argparse
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import new_majiq as nm
import new_majiq.constants as constants
from new_majiq._run._majiq_args import check_nonnegative_factory
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

# type aliases
SJGroupT = Sequence[Path]
SJGroupsT = Dict[str, SJGroupT]


# what logic do we follow?
class BuildType(Enum):
    BUILD_ALL = "build_all"
    BUILD_KNOWN = "build_known"
    SIMPLIFY_ONLY = "simplify_only"


class IntronsType(Enum):
    NO_INTRONS = "no_introns"
    ANNOTATED_INTRONS = "annotated_introns"
    ALL_INTRONS = "all_introns"


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
        "--reset-simplified",
        action="store_true",
        dest="reset_simplify",
        default=True,
        help="If simplifying, reset all introns/junctions to simplified,"
        " i.e. start simplification from scratch"
        " (default: reset_simplify=%(default)s)",
    )
    reset_ex.add_argument(
        "--update-simplified",
        action="store_false",
        dest="reset_simplify",
        default=True,
        help="Continue from existing simplified splicegraph,"
        " i.e. new denovo connections are simplified, existing connections will"
        " not be reset to simplified (default: reset_simplify=%(default)s)",
    )
    return


def build_threshold_args(parser: argparse.ArgumentParser) -> None:
    """arguments for build threshold parameters"""
    thresholds = parser.add_argument_group("Build filters")
    # min-experiments
    thresholds.add_argument(
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MINEXPERIMENTS,
        metavar="x",
        help="Threshold per experiments group. If < 1, the fraction of"
        " experiments in a group that must pass individual filters for a"
        " feature to be accepted. If greater, an absolute number."
        " (default: %(default)s)",
    )
    # per-experiment thresholds
    thresholds.add_argument(
        "--minreads",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINREADS,
        metavar="x",
        help="Minimum readrate to pass an annotated junction (default: %(default)s)",
    )
    thresholds.add_argument(
        "--mindenovo",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINDENOVO,
        metavar="x",
        help="Minimum readrate to pass a denovo junction or intron (default: %(default)s)",
    )
    thresholds.add_argument(
        "--minpos",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BUILD_MINPOS,
        metavar="x",
        help="Minimum number of nonzero positions to pass a junction."
        " This is scaled for introns with some minimum coverage per bin to"
        " account for length dependence."
        " (default: %(default)s).",
    )
    thresholds.add_argument(
        "--max-pctbins",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MAX_PCTBINS,
        metavar="x",
        help="Maximum fraction of intron bins on which to require coverage"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        "--junction-acceptance-probability",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
        help="Set length-dependent minbins intron thresholds by considering"
        " per-position Poisson readrate for which junctions are accepted with"
        " this probability (default: %(default)s)",
    )
    thresholds.add_argument(
        "--intron-acceptance-probability",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_BUILD_MATCH_INTRON_PROBABILITY,
        help="Set length-dependent minbins intron thresholds by picking least"
        " strict thresholds for which per-position readrate determined by"
        " --junction-acceptance-probability gives has acceptance probability"
        " less than this probability (default: %(default)s)",
    )
    return


def ir_filtering_args(parser: argparse.ArgumentParser) -> None:
    """arguments for intron filtering"""
    introns = parser.add_argument_group("Intron filtering")
    introns_ex = introns.add_mutually_exclusive_group()
    introns_ex.add_argument(
        "--all-introns",
        dest="introns",
        default=IntronsType.ALL_INTRONS,
        action="store_const",
        const=IntronsType.ALL_INTRONS,
        help="Keep all annotated introns and denovo introns passing build filters"
        " (default: %(default)s)",
    )
    introns_ex.add_argument(
        "--annotated-introns",
        dest="introns",
        default=IntronsType.ALL_INTRONS,
        action="store_const",
        const=IntronsType.ANNOTATED_INTRONS,
        help="Keep all annotated introns only (default: %(default)s)",
    )
    introns_ex.add_argument(
        "--no-introns",
        dest="introns",
        default=IntronsType.ALL_INTRONS,
        action="store_const",
        const=IntronsType.NO_INTRONS,
        help="Drop/ignore all introns (default: %(default)s)",
    )
    return


def build_type_args(parser: argparse.ArgumentParser) -> None:
    """arguments for build type"""
    build = parser.add_argument_group("Build type")
    build_ex = build.add_mutually_exclusive_group()
    build_ex.add_argument(
        "--build-all",
        dest="build",
        default=BuildType.BUILD_ALL,
        action="store_const",
        const=BuildType.BUILD_ALL,
        help="Process experiments to find new junctions and update existing junctions"
        " (default: %(default)s)",
    )
    build_ex.add_argument(
        "--build-known",
        dest="build",
        default=BuildType.BUILD_ALL,
        action="store_const",
        const=BuildType.BUILD_KNOWN,
        help="Process experiments to update known junctions"
        " (note that denovo/annotated intron processing specified separately)"
        " (default: %(default)s)",
    )
    build_ex.add_argument(
        "--simplify-only",
        dest="build",
        default=BuildType.BUILD_ALL,
        action="store_const",
        const=BuildType.SIMPLIFY_ONLY,
        help="Only perform simplification (default: %(default)s)",
    )
    return


def enable_simplify_args(parser: argparse.ArgumentParser) -> None:
    simplify = parser.add_argument_group("Simplification")
    simplify_ex = simplify.add_mutually_exclusive_group()
    simplify_ex.add_argument(
        "--simplify",
        dest="simplify",
        action="store_true",
        default=None,
        help="(Un)simplify splicegraph using evidence from input experiments"
        " (default: do not simplify unless build type set to --simplify-only)",
    )
    simplify_ex.add_argument(
        "--no-simplify",
        dest="simplify",
        action="store_false",
        default=None,
        help="Explicitly request to not perform simplification."
        " Will raise error if --simplify-only is set."
        " (default: do not simplify unless build type set to --simplify-only)",
    )
    return


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument("base_sg", type=Path, help="Path to base splicegraph")
    parser.add_argument("out_sg", type=Path, help="Path for output splicegraph")
    experiments = parser.add_argument_group("Input experiments (one required)")
    experiments_ex = experiments.add_mutually_exclusive_group(required=True)
    experiments_ex.add_argument(
        "--groups-tsv",
        type=Path,
        metavar="TSV",
        help="Specify experiments from multiple build groups using TSV file."
        " Required columns 'group' and 'sj'."
        " One row per unique experiment."
        " `group` indicates group experiment belongs to, `sj`"
        " the path to the experiment's SJ file (from `new-majiq sj`)",
    )
    experiments_ex.add_argument(
        "--sjs",
        type=Path,
        nargs="+",
        help="Specify experiments from a single build group directly as"
        " the unique paths to the group experiments' SJ files"
        " (use `new-majiq combine` to merge independent build groups)",
    )

    build_type_args(parser)
    ir_filtering_args(parser)
    enable_simplify_args(parser)

    reset_simplified_args(parser)
    build_threshold_args(parser)
    simplifier_threshold_args(parser, prefix="simplify-")
    return


def get_grouped_experiments(groups_tsv: Path) -> SJGroupsT:
    """Get grouped experiments from table

    Get grouped experiments from table. Verify that there are no duplicate
    experiments (insofar as their paths) and that the paths exist.
    """
    import pandas as pd

    df = pd.read_csv(groups_tsv, sep="\t", usecols=["group", "sj"])
    if len(df.sj) == 0:
        raise ValueError(f"No input experiments specified in {groups_tsv}")
    # determine if any duplicated experiments
    duplicated_mask = df.sj.duplicated()
    if duplicated_mask.any():
        duplicated_sj = set(df.sj[duplicated_mask])
        raise ValueError(f"Requested build with repeated experiments {duplicated_sj}")
    # verify that all paths exist
    for sj_path in df.sj:
        if not Path(sj_path).exists():
            raise ValueError(f"Unable to find input experiment {sj_path}")
    return {
        group: sorted(Path(x) for x in group_sjs)
        for group, group_sjs in df.groupby("group")["sj"]
    }


def do_build(
    sg: nm.SpliceGraph,
    experiments: SJGroupsT,
    experiment_thresholds: nm.ExperimentThresholds,
    min_experiments: float,
    process_denovo_junctions: bool,
    introns_type: IntronsType,
    denovo_simplified: bool,
) -> nm.SpliceGraph:
    log = get_logger()

    log.info("Updating known and identifying denovo junctions")
    junction_builder = sg.junctions.builder()
    for group_ndx, (group, group_sjs) in enumerate(experiments.items()):
        log.info(
            f"Processing junctions from group {group} ({1 + group_ndx} / {len(experiments)})"
        )
        build_group = sg.junctions.build_group(sg.exons)
        for sj_ndx, sj in enumerate(group_sjs):
            log.info(
                f"Processing junctions from {Path(sj).resolve()}"
                f" (experiment {1 + sj_ndx} / {len(group_sjs)} in group)"
            )
            build_group.add_experiment(
                nm.SJJunctionsBins.from_zarr(sj),
                thresholds=experiment_thresholds,
                add_denovo=process_denovo_junctions,
            )
        junction_builder.add_group(build_group, min_experiments)
    del build_group
    log.info("Consolidating passed junctions from input experiments")
    updated_junctions = junction_builder.get_passed(denovo_simplified)
    del junction_builder

    log.info("Inferring denovo exons and updated exon boundaries")
    updated_exons = sg.exons.infer_with_junctions(updated_junctions)

    updated_introns: nm.GeneIntrons
    if introns_type == IntronsType.NO_INTRONS:
        log.info("Dropping all introns")
        updated_introns = updated_exons.empty_introns()
    else:
        log.info("Determining potential gene introns using updated exons")
        potential_introns = updated_exons.potential_introns(denovo_simplified)
        potential_introns.update_flags_from(sg.introns)
        log.info("Identiying new passed introns")
        intron_group = potential_introns.build_group()  # intron groups done in place
        for group_ndx, (group, group_sjs) in enumerate(experiments.items()):
            log.info(
                f"Processing introns from group {group}"
                f" ({1 + group_ndx} / {len(experiments)})"
            )
            for sj_ndx, sj in enumerate(group_sjs):
                log.info(
                    f"Processing introns from {Path(sj).resolve()}"
                    f" (experiment {1 + sj_ndx} / {len(group_sjs)} in group)"
                )
                intron_group.add_experiment(
                    nm.SJIntronsBins.from_zarr(sj),
                    thresholds=experiment_thresholds,
                )
            intron_group.update_introns(min_experiments)
        del intron_group
        log.info("Filtering potential introns to those passing thresholds")
        updated_introns = potential_introns.filter_passed(
            keep_annotated=True,
            discard_denovo=introns_type != IntronsType.ALL_INTRONS,
        )
        del potential_introns

    log.info("constructing splicegraph with updated exons and connections")
    return sg.with_updated_exon_connections(
        nm.ExonConnections.create_connecting(
            updated_exons, updated_introns, updated_junctions
        )
    )


def do_simplify(
    sg: nm.SpliceGraph,
    experiments: SJGroupsT,
    reset_simplify: bool,
    simplify_minpsi: float,
    simplify_minreads_annotated: float,
    simplify_minreads_denovo: float,
    simplify_minreads_ir: float,
    simplify_min_experiments: float,
) -> nm.SpliceGraph:
    log = get_logger()

    if reset_simplify:
        log.info("Setting all introns and junctions to simplified")
        sg.introns._simplify_all()
        sg.junctions._simplify_all()

    log.info("Identifying introns and junctions to unsimplify")
    simplifier_group = sg.exon_connections.simplifier()
    for group_ndx, (group, group_sjs) in enumerate(experiments.items()):
        log.info(
            f"Processing coverage from group {group} ({1 + group_ndx} / {len(experiments)})"
        )
        for sj_ndx, sj in enumerate(group_sjs):
            log.info(
                f"Processing coverage from {Path(sj).resolve()}"
                f" (experiment {1 + sj_ndx} / {len(group_sjs)} in group)"
            )
            simplifier_group.add_experiment(
                nm.SJExperiment.from_zarr(sj),
                min_psi=simplify_minpsi,
                minreads_annotated=simplify_minreads_annotated,
                minreads_denovo=simplify_minreads_denovo,
                minreads_introns=simplify_minreads_ir,
            )
        simplifier_group.update_connections(simplify_min_experiments)
    del simplifier_group

    return sg


def run(args: argparse.Namespace) -> None:
    simplify: bool  # will we be simplifying in the end?
    if args.simplify is None:
        simplify = args.build == BuildType.SIMPLIFY_ONLY
    else:
        if not args.simplify and args.build == BuildType.SIMPLIFY_ONLY:
            raise argparse.ArgumentError(
                None, "--no-simplify and --simplify-only are incompatible"
            )
        simplify = args.simplify

    if not args.base_sg.exists():
        raise ValueError(f"Unable to find base splicegraph at {args.base_sg}")
    if args.out_sg.exists():
        raise ValueError(f"Output path {args.out_sg} already exists")

    log = get_logger()
    experiment_thresholds = nm.ExperimentThresholds(
        minreads=args.minreads,
        mindenovo=args.mindenovo,
        minpos=args.minpos,
        max_pctbins=args.max_pctbins,
        junction_acceptance_probability=args.junction_acceptance_probability,
        intron_acceptance_probability=args.intron_acceptance_probability,
    )
    if args.groups_tsv:
        if not args.groups_tsv.exists():
            raise ValueError(
                f"Unable to find group definitions at {args.groups_tsv.resolve()}"
            )
        log.info("Loading experiment groups")
        experiments = get_grouped_experiments(args.groups_tsv)
    else:
        if missing := sorted(set(x for x in args.sjs if not x.exists())):
            raise ValueError(f"Unable to find all input SJ files ({missing =})")
        if len(unique := set(args.sjs)) != len(args.sjs):
            # get non-unique sjs to report error
            non_unique = sorted(args.sjs)
            for x in unique:
                non_unique.remove(x)  # removes first occurence
            raise ValueError(f"Non-unique input SJ files ({non_unique = })")
        experiments = {"": args.sjs}
    log.info(f"Loading base splicegraph from {args.base_sg.resolve()}")
    sg = nm.SpliceGraph.from_zarr(args.base_sg)

    # perform build?
    if args.build != BuildType.SIMPLIFY_ONLY:
        sg = do_build(
            sg=sg,
            experiments=experiments,
            experiment_thresholds=experiment_thresholds,
            min_experiments=args.min_experiments,
            # we are procesing denovo junction only if build_all specified
            process_denovo_junctions=args.build == BuildType.BUILD_ALL,
            introns_type=args.introns,
            # if simplifying or if --update-simplified (i.e. not --reset-simplified)
            denovo_simplified=simplify or not args.reset_simplify,
        )
    else:
        log.info(f"Filtering introns {args.introns = }")
        sg = nm.SpliceGraph.from_components(
            contigs=sg.contigs,
            genes=sg.genes,
            exons=sg.exons,
            junctions=sg.junctions,
            introns=(
                sg.introns.filter_passed(keep_annotated=True, discard_denovo=False)
                if args.introns == IntronsType.ALL_INTRONS
                else sg.introns.filter_passed(keep_annotated=True, discard_denovo=True)
                if args.introns == IntronsType.ANNOTATED_INTRONS
                else sg.exons.empty_introns()
            ),
        )

    # perform simplification?
    if simplify:
        sg = do_simplify(
            sg=sg,
            experiments=experiments,
            reset_simplify=args.reset_simplify,
            simplify_minpsi=args.simplify_minpsi,
            simplify_minreads_annotated=args.simplify_minreads_annotated,
            simplify_minreads_denovo=args.simplify_minreads_denovo,
            simplify_minreads_ir=args.simplify_minreads_ir,
            simplify_min_experiments=args.simplify_min_experiments,
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
