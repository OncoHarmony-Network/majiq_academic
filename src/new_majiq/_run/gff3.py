"""
gff3.py

Translate input GFF3 file to base annotated splicegraph file

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import ExistingResolvedPath, NewResolvedPath
from new_majiq._run._run import GenericSubcommand
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
    # annotate ir options
    introns_ex = parser.add_mutually_exclusive_group()
    introns_ex.add_argument(
        "--ignore_ir",
        action="store_false",
        dest="process_ir",
        default=nm.constants.DEFAULT_BUILD_PROCESS_IR,
        help="Ignore annotated introns when processing annotations"
        " (default process_ir=%(default)s)",
    )
    introns_ex.add_argument(
        "--process_ir",
        action="store_true",
        dest="process_ir",
        default=nm.constants.DEFAULT_BUILD_PROCESS_IR,
        help="Ensure annotated introns are processed (default process_ir=%(default)s)",
    )
    # TODO(jaicher): enable configuration of GFF3 parsing
    return


def run(args: argparse.Namespace) -> None:
    log = get_logger()
    log.info(f"Processing GFF3 {args.gff3} to create annotated splicegraph")
    sg = nm.SpliceGraph.from_gff3(
        args.gff3,
        process_ir=args.process_ir,
        # TODO(jaicher): enable configuration of GFF3 parsing
        gff3_types=nm.constants.DEFAULT_BUILD_GFF3TYPES,
    )
    # TODO(jaicher): enable capture of skipped features to print to logger
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
