"""
combine.py

Combine input splicegraphs

Author: Joseph K Aicher
"""

import argparse

from new_majiq._run.build import ir_filtering_args
from pathlib import Path
from new_majiq._run._run import GenericSubcommand
from typing import (
    List,
    Optional,
)


DESCRIPTION = "Combine input splicegraphs into single splicegraph"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "out_sg",
        type=Path,
        help="Path for output splicegraph",
    )
    input_args = parser.add_argument_group("Input splicegraphs (need at least one)")
    input_args.add_argument(
        "--make-annotated",
        type=Path,
        nargs="+",
        default=[],
        help="Input splicegraphs for which all junctions will be marked as"
        " annotated (i.e. not denovo). This helps highlight denovo junctions"
        " that were unique to non-base splicegraphs. Note that introns remain"
        " denovo in order to accommodate assignment of intronic coverage.",
    )
    input_args.add_argument(
        "--keep-denovo",
        type=Path,
        nargs="+",
        default=list(),
        help="Input splicegraphs for which junctions will remain marked as"
        " denovo (unless also found in inputs from --make-annotated)",
    )

    ir_filtering_args(parser)
    return


def run(args: argparse.Namespace) -> None:
    if args.out_sg.exists():
        raise ValueError(
            f"Output splicegraph {args.output_sg.resolve()} already exists"
        )
    all_inputs = args.make_annotated + args.keep_denovo
    if not all_inputs:
        raise ValueError("No input splicegraphs were provided to combine")
    if (missing := sorted(set(x for x in all_inputs if not x.exists()))) :
        raise ValueError(f"Unable to find all input splicegraphs ({missing = })")

    import new_majiq as nm
    from new_majiq.Genes import Genes
    from new_majiq.GeneJunctions import GeneJunctions
    from new_majiq.Exons import Exons
    from new_majiq.GeneIntrons import GeneIntrons
    from new_majiq.logger import get_logger

    log = get_logger()

    log.info("Combining input junctions")
    df_junctions = GeneJunctions.combine_datasets(
        [
            GeneJunctions.load_dataset(p).assign_coords(denovo=False)
            for p in args.make_annotated
        ]
        + [
            GeneJunctions.load_dataset(p)
            for p in args.keep_denovo
        ]
    )
    log.info("Constructing GeneJunctions object for combined set of junctions")
    genes = Genes.from_zarr(all_inputs[0])
    junctions = GeneJunctions.from_dataset_and_genes(df_junctions, genes)
    log.info("Creating updated combined exon definitions")
    exons = Exons.from_zarr(all_inputs[0], genes).infer_with_junctions(junctions)
    log.info("Defining new potential introns between exons")
    potential_introns = exons.potential_introns(True)  # all start simplified
    log.info("Updating intron flags using input splicegraphs")
    for p in all_inputs:
        potential_introns.update_flags_from(GeneIntrons.from_zarr(p, genes))
    log.info("Filtering introns to those passing thresholds")
    introns = potential_introns.filter_passed(
        keep_annotated=args.keep_annotated_ir,
        discard_denovo=not args.process_denovo_introns,
    )
    del potential_introns
    log.info("Creating final combined splicegraph")
    sg = nm.SpliceGraph.from_components(genes.contigs, genes, exons, junctions, introns)
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
