"""
build_group.py

Update splicegraph with a single build group

Author: Joseph K Aicher
"""

import argparse

from pathlib import Path
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build import (
    build_threshold_args,
    denovo_junctions_args,
    ir_filtering_args,
    denovo_simplified_args,
)

from typing import (
    List,
    Optional,
)


DESCRIPTION = "Update splicegraph with a single group of experiments"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument("base_sg", type=Path, help="Path to base splicegraph")
    parser.add_argument("out_sg", type=Path, help="Path for output splicegraph")
    parser.add_argument(
        "sj", type=Path, nargs="+", help="Paths to input experiments as SJ files"
    )
    build_threshold_args(parser)
    denovo_junctions_args(parser)
    ir_filtering_args(parser)
    denovo_simplified_args(parser)
    return


def run(args: argparse.Namespace) -> None:
    if not args.base_sg.exists():
        raise ValueError(f"Unable to find base splicegraph at {args.base_sg}")
    if (missing := sorted(set(x for x in args.sj if not x.exists()))):
        raise ValueError(f"Unable to find all input SJ files ({missing =})")
    if args.out_sg.exists():
        raise ValueError(f"Output path {args.out_sg} already exists")

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

    log.info(f"Loading base splicegraph from {args.base_sg.resolve()}")
    sg = nm.SpliceGraph.from_zarr(args.base_sg)

    log.info("Updating known and identifying denovo junctions")
    junction_builder = sg.junctions.builder()
    build_group = sg.junctions.build_group(sg.exons)
    for ndx, sj in enumerate(args.sj):
        log.info(
            f"Processing junctions from {sj.resolve()}"
            f" (experiment {1 + ndx} / {len(args.sj)})"
        )
        build_group.add_experiment(
            nm.SJJunctionsBins.from_zarr(sj),
            thresholds=experiment_thresholds,
            add_denovo=args.process_denovo_junctions,
        )
    log.info("Consolidating passed junctions from the build group")
    junction_builder.add_group(build_group, args.min_experiments)
    updated_junctions = junction_builder.get_passed(args.denovo_simplified)
    del build_group, junction_builder

    log.info("Inferring denovo exons and updated exon boundaries")
    updated_exons = sg.exons.infer_with_junctions(updated_junctions)

    log.info("Determining potential gene introns using updated exons")
    potential_introns = updated_exons.potential_introns(args.denovo_simplified)
    potential_introns.update_flags_from(sg.introns)
    log.info("Identiying new passed introns")
    intron_group = potential_introns.build_group()  # intron groups done in place
    for ndx, sj in enumerate(args.sj):
        log.info(
            f"Processing introns from {sj.resolve()}"
            f" (experiment {1 + ndx} / {len(args.sj)})"
        )
        intron_group.add_experiment(
            nm.SJIntronsBins.from_zarr(sj),
            thresholds=experiment_thresholds,
        )
    intron_group.update_introns(args.min_experiments)
    del intron_group
    log.info("Filtering potential introns to those passing thresholds")
    updated_introns = potential_introns.filter_passed(
        keep_annotated=args.keep_annotated_ir,
        discard_denovo=not args.process_denovo_introns,
    )
    del potential_introns

    log.info("constructing splicegraph with updated exons and connections")
    sg = sg.with_updated_exon_connections(
        ExonConnections.create_connecting(
            updated_exons, updated_introns, updated_junctions
        )
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
