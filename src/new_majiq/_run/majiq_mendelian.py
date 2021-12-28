#!/usr/bin/env python
"""
majiq_mendelian.py

Provides entry-point to majiq scripts for assessing patients with suspected
Mendelian disorders

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict

from new_majiq._run._run import GenericSubcommand
from new_majiq._run.psi_controls import subcommand as psi_controls
from new_majiq._run.psi_outliers import subcommand as psi_outliers
from new_majiq._run.psi_outliers_summary import subcommand as psi_outliers_summary

SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "controls": psi_controls,
    "outliers": psi_outliers,
    "outliers-summary": psi_outliers_summary,
}


def main() -> None:
    """Entry-point into multiple tools using subcommands"""
    # build parser
    parser = argparse.ArgumentParser(
        description="Tools to summarize PSI in controls and identify outliers in cases"
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
