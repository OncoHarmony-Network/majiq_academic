"""
_run.py

Shared types, functions for majiq-tools subcommands

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Callable
from new_majiq.logger import setup_logger, get_logger


class GenericSubcommand(object):
    """ Wrap add_args, run with shared setup (especially for logging)
    """

    def __init__(
        self,
        DESCRIPTION: str,
        add_args: Callable[[argparse.ArgumentParser], None],
        run: Callable[[argparse.Namespace], None],
    ) -> None:
        self._description = DESCRIPTION
        self._add_args = add_args
        self._run = run
        return

    @property
    def DESCRIPTION(self) -> str:
        return self._description

    def add_args(self, parser: argparse.ArgumentParser) -> None:
        """ Add arguments to provided argument parser
        """
        self._add_args(parser)
        # add parameters for logging
        parser.add_argument(
            "--logger",
            type=str,
            default=None,
            help="Redirect logging to specified file (default: stderr)",
        )
        parser.add_argument(
            "--silent", action="store_true", default=False, help="Silence the logger"
        )
        parser.add_argument(
            "--debug",
            action="store_true",
            default=False,
            help="Enable detailed logging for debugging purposes",
        )
        return

    def run(self, args: argparse.Namespace) -> None:
        """ Run subcommand with parsed arguments
        """
        # set up logging
        setup_logger(logfile=args.logger, silent=args.silent, debug=args.debug)
        log = get_logger()
        # print information about the run
        from new_majiq.version import version
        log.info(f"new-majiq v{version}")
        log.info(f"Command: {' '.join(sys.argv)}")
        log.info(f"Arguments:\n{args}")
        # run subcommand
        try:
            self._run(args)
        except Exception:
            # log exception and re-raise
            log.exception("Exiting due to exception:")
            sys.exit(-1)
            return
        # when done running, note that it was successsfully completed!
        log.info("Finished successfully!")
        return
