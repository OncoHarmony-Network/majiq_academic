"""
gff3.py

Translate input GFF3 file to base annotated splicegraph file

Author: Joseph K Aicher
"""

import argparse

import new_majiq.constants as constants

from pathlib import Path
from new_majiq._run._run import GenericSubcommand
from typing import (
    List,
    Optional,
)


DESCRIPTION = "Translate input GFF3 file to base splicegraph file"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "gff3", type=Path, help="Path to GFF3 file (uncompressed or gzipped) to convert"
    )
    parser.add_argument(
        "splicegraph",
        type=Path,
        help="Path to save resulting splicegraph. Fails if path already exists.",
    )
    # annotate ir options
    introns_ex = parser.add_mutually_exclusive_group()
    introns_ex.add_argument(
        "--ignore_ir",
        action="store_false",
        dest="process_ir",
        default=constants.DEFAULT_BUILD_PROCESS_IR,
        help="Ignore annotated introns when processing annotations"
        " (default process_ir=%(default)s)",
    )
    introns_ex.add_argument(
        "--process_ir",
        action="store_true",
        dest="process_ir",
        default=constants.DEFAULT_BUILD_PROCESS_IR,
        help="Ensure annotated introns are processed (default process_ir=%(default)s)",
    )
    # TODO(jaicher): enable configuration of GFF3 parsing
    return


def run(args: argparse.Namespace) -> None:
    if not args.gff3.exists():
        raise ValueError(f"Was unable to find input GFF3 file at {args.gff3}")
    if args.splicegraph.exists():
        raise ValueError(
            f"Path for output splicegraph ({args.splicegraph}) already exists"
        )
    from new_majiq import SpliceGraph
    from new_majiq.logger import get_logger

    log = get_logger()
    log.info(f"Processing GFF3 {args.gff3.resolve()} to create annotated splicegraph")
    sg = SpliceGraph.from_gff3(
        args.gff3,
        process_ir=args.process_ir,
        # TODO(jaicher): enable configuration of GFF3 parsing
        gff3_types=constants.DEFAULT_BUILD_GFF3TYPES,
    )
    # TODO(jaicher): enable capture of skipped features to print to logger
    log.info(f"Saving annotated splicegraph to {args.splicegraph.resolve()}")
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
