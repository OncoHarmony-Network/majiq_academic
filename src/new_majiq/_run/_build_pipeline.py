"""
_build_pipeline.py

Components of build pipeline for use in constructing scripts/subcommands

Author: Joseph K Aicher
"""

import new_majiq as nm

from pathlib import Path
from new_majiq.ExonConnections import ExonConnections
from new_majiq.logger import get_logger
from typing import (
    Dict,
    Sequence,
    Union,
)
from enum import Enum


# type aliases
SJGroupT = Sequence[Path]
SJGroupsT = Dict[str, SJGroupT]


class BuildType(Enum):
    BUILD_ALL = "build_all"
    BUILD_KNOWN = "build_known"
    SIMPLIFY_ONLY = "simplify_only"


class IntronsType(Enum):
    NO_INTRONS = "no_introns"
    ANNOTATED_INTRONS = "annotated_introns"
    ALL_INTRONS = "all_introns"


def get_grouped_experiments(groups_tsv: Path) -> SJGroupsT:
    """Get grouped experiments from table

    Get grouped experiments from table. Verify that there are no duplicate
    experiments (insofar as their paths) and that the paths exist.
    """
    import pandas as pd

    df = pd.read_csv(groups_tsv, sep="\t", usecols=["group", "sj"])
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


def build(
    sg: Union[Path, nm.SpliceGraph],
    experiments: SJGroupsT,
    experiment_thresholds: nm.ExperimentThresholds,
    min_experiments: float,
    process_denovo_junctions: bool,
    introns_type: IntronsType,
    denovo_simplified: bool,
) -> nm.SpliceGraph:
    log = get_logger()
    if not isinstance(sg, nm.SpliceGraph):
        log.info(f"Loading base splicegraph from {sg.resolve()}")
        sg = nm.SpliceGraph.from_zarr(sg)

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
        ExonConnections.create_connecting(
            updated_exons, updated_introns, updated_junctions
        )
    )


def simplify(
    sg: Union[Path, nm.SpliceGraph],
    experiments: SJGroupsT,
    reset_simplify: bool,
    simplify_minpsi: float,
    simplify_minreads_annotated: float,
    simplify_minreads_denovo: float,
    simplify_minreads_ir: float,
    simplify_min_experiments: float,
) -> nm.SpliceGraph:
    log = get_logger()
    if not isinstance(sg, nm.SpliceGraph):
        log.info(f"Loading base splicegraph from {sg.resolve()}")
        sg = nm.SpliceGraph.from_zarr(sg)

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
                nm.SpliceGraphReads.from_connections_and_sj(
                    sg.introns,
                    sg.junctions,
                    nm.SJIntronsBins.from_zarr(sj),
                    nm.SJJunctionsBins.from_zarr(sj),
                ),
                min_psi=simplify_minpsi,
                minreads_annotated=simplify_minreads_annotated,
                minreads_denovo=simplify_minreads_denovo,
                minreads_introns=simplify_minreads_ir,
            )
        simplifier_group.update_connections(simplify_min_experiments)
    del simplifier_group

    return sg
