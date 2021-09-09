#!/usr/bin/env python
"""
cite.py
"""

import argparse
from typing import List, Optional

from new_majiq._run._run import GenericSubcommand

DESCRIPTION = "Information about how to cite MAJIQ"
MESSAGE = r"""
 /$$      /$$  /$$$$$$     /$$$$$ /$$$$$$  /$$$$$$
| $$$    /$$$ /$$__  $$   |__  $$|_  $$_/ /$$__  $$
| $$$$  /$$$$| $$  \ $$      | $$  | $$  | $$  \ $$
| $$ $$/$$ $$| $$$$$$$$      | $$  | $$  | $$  | $$
| $$  $$$| $$| $$__  $$ /$$  | $$  | $$  | $$  | $$
| $$\  $ | $$| $$  | $$| $$  | $$  | $$  | $$/$$ $$
| $$ \/  | $$| $$  | $$|  $$$$$$/ /$$$$$$|  $$$$$$/
|__/     |__/|__/  |__/ \______/ |______/ \____ $$$
                                               \__/

Please cite us! Recommended DOIs:
+ (TODO)
"""


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser -- but we have none to add"""
    return


def run(args: argparse.Namespace) -> None:
    print(MESSAGE)
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
