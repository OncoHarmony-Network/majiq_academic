"""
combine.py

Combine input splicegraphs

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build_args import IntronsType, ir_filtering_args
from new_majiq.logger import get_logger

DESCRIPTION = "Combine input splicegraphs into single splicegraph"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "out_sg",
        type=NewResolvedPath,
        help="Path for output splicegraph",
    )
    StoreSGPaths = StoreRequiredUniqueActionFactory()
    input_args = parser.add_argument_group(
        "required splicegraph arguments (need at least one)"
    )
    input_args.add_argument(
        "--make-annotated",
        metavar="SG",
        type=ExistingResolvedPath,
        action=StoreSGPaths,
        nargs="+",
        default=list(),
        help="Input splicegraphs for which all junctions will be marked as"
        " annotated (i.e. not denovo). This helps highlight denovo junctions"
        " that were unique to non-base splicegraphs. Note that introns remain"
        " denovo in order to accommodate assignment of intronic coverage.",
    )
    input_args.add_argument(
        "--keep-denovo",
        metavar="SG",
        type=ExistingResolvedPath,
        action=StoreSGPaths,
        nargs="+",
        default=list(),
        help="Input splicegraphs for which junctions will remain marked as"
        " denovo (unless also found in inputs from --make-annotated)",
    )

    ir_filtering_args(parser)
    return


def run(args: argparse.Namespace) -> None:
    all_inputs = args.make_annotated + args.keep_denovo
    if not all_inputs:
        raise ValueError("No input splicegraphs were provided to combine")

    log = get_logger()

    log.info("Combining input junctions")
    genes = nm.Genes.from_zarr(all_inputs[0])
    accumulator = nm.GeneJunctionsAccumulator(genes)
    for p in args.make_annotated:
        accumulator.add(nm.GeneJunctions.from_zarr(p, genes=genes), make_annotated=True)
    for p in args.keep_denovo:
        accumulator.add(
            nm.GeneJunctions.from_zarr(p, genes=genes), make_annotated=False
        )
    junctions = accumulator.accumulated()
    log.info("Creating updated combined exon definitions")
    exons = nm.Exons.from_zarr(all_inputs[0], genes).infer_with_junctions(junctions)
    introns: nm.GeneIntrons
    if args.introns == IntronsType.NO_INTRONS:
        log.info("Creating matching empty introns")
        introns = exons.empty_introns()
    else:
        log.info("Defining new potential introns between exons")
        introns = exons.potential_introns(True)  # all start simplified
        potential_introns = exons.get_annotated().potential_introns(True)
        log.info("Updating intron flags using input splicegraphs")
        # propagate intron flags within annotated exon boundaries
        for p in all_inputs:
            potential_introns.update_flags_from(nm.GeneIntrons.from_zarr(p, genes))
        # propagate back to introns between updated exon boundaries
        potential_introns = potential_introns.filter_passed(
            keep_annotated=True,
            discard_denovo=args.introns == IntronsType.ANNOTATED_INTRONS,
        )
        introns.update_flags_from(potential_introns)
        del potential_introns  # clear from memory, no longer necessary
        log.info("Filtering introns to those passing thresholds")
        introns = introns.filter_passed(
            keep_annotated=True,
            discard_denovo=args.introns == IntronsType.ANNOTATED_INTRONS,
        )
    log.info("Creating final combined splicegraph")
    sg = nm.SpliceGraph.from_components(genes.contigs, genes, exons, junctions, introns)
    log.info(f"Saving updated splicegraph to {args.out_sg}")
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
