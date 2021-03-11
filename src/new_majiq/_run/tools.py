#!/usr/bin/env python
"""
tools.py

Provides entry-point to different majiq-tools scripts

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Optional, List, Dict
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.cite import subcommand as cite
from new_majiq._run.gff3 import subcommand as gff3


SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "gff3": gff3,
    "cite": cite,
}


def main(sys_args: Optional[List[str]] = None) -> None:
    """ Entry-point into multiple tools using subcommands
    """
    # build parser
    parser = argparse.ArgumentParser(
        description="Tools to detect, quantify, and analyze RNA splicing"
    )
    # add subparsers
    subparsers = parser.add_subparsers(required=True, help="")
    for src_name, src_module in SUBPARSER_SOURCES.items():
        src_parser = subparsers.add_parser(
            src_name, help=src_module.DESCRIPTION, description=src_module.DESCRIPTION,
        )
        src_parser.set_defaults(func=src_module.run)
        src_module.add_args(src_parser)

    # check length of input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # parse arguments now
    args = parser.parse_args(sys_args)
    args.func(args)


if __name__ == "__main__":
    main()
