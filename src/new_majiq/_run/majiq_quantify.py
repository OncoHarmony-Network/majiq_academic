#!/usr/bin/env python
"""
majiq_quantify.py

Provides entry-point to different majiq quantifiers

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict

from new_majiq import __version__
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.deltapsi import subcommand as deltapsi
from new_majiq._run.heterogen import subcommand as heterogen
from new_majiq._run.quantify import subcommand as quantify

SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "psi": quantify,
    "deltapsi": deltapsi,
    "heterogen": heterogen,
}


def main() -> None:
    """Entry-point into multiple tools using subcommands"""
    # build parser
    parser = argparse.ArgumentParser(
        description="Tools to quantify RNA splicing within or between groups"
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    # add subparsers
    subparsers = parser.add_subparsers(required=True, help="")
    for src_name, src_module in SUBPARSER_SOURCES.items():
        src_parser = subparsers.add_parser(
            src_name,
            help=src_module.DESCRIPTION,
            description=src_module.DESCRIPTION,
        )
        src_parser.set_defaults(func=src_module.run)
        src_module.add_args(src_parser)

    # check length of input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # parse arguments now
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
