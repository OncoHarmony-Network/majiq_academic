"""
gff3.py

Translate input GFF3 file to base annotated splicegraph file

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

from new_majiq._run._majiq_args import ExistingResolvedPath, NewResolvedPath
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build_args import gff3_parsing_args, gff3_parsing_pipeline
from new_majiq.logger import get_logger

DESCRIPTION = "Translate input GFF3 file to base splicegraph file"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "gff3",
        type=ExistingResolvedPath,
        help="Path to GFF3 file (uncompressed or gzipped) to convert",
    )
    parser.add_argument(
        "splicegraph",
        type=NewResolvedPath,
        help="Path to save resulting splicegraph. Fails if path already exists.",
    )
    # enable configuration of GFF3 parsing
    gff3_parsing_args(parser)
    return


def run(args: argparse.Namespace) -> None:
    log = get_logger()
    sg = gff3_parsing_pipeline(
        args.gff3,
        features_default=args.gff3_features_default,
        features_tsv=args.gff3_features_tsv,
        types_genes=args.gff3_types_genes,
        types_transcripts=args.gff3_types_transcripts,
        types_exons=args.gff3_types_exons,
        types_silent=args.gff3_types_silent,
        types_hard_skip=args.gff3_types_hard_skip,
    )
    log.info(f"Saving annotated splicegraph to {args.splicegraph}")
    sg.to_zarr(args.splicegraph)
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
