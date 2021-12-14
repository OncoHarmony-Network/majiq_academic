"""
_run.py

Shared types, functions for majiq-tools subcommands

Author: Joseph K Aicher
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Callable

from dask.distributed import Client

from new_majiq.logger import get_logger, setup_logger


class GenericSubcommand(object):
    """Wrap add_args, run with shared setup (especially for logging)"""

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
        """Add arguments to provided argument parser"""
        self._add_args(parser)
        # add parameters for logging
        log_args = parser.add_argument_group("logging arguments")
        log_args.add_argument(
            "--logger",
            type=str,
            default=None,
            help="Redirect logging to specified file (default: stderr)",
        )
        log_args.add_argument(
            "--silent", action="store_true", default=False, help="Silence the logger"
        )
        log_args.add_argument(
            "--debug",
            action="store_true",
            default=False,
            help="Enable detailed logging for debugging purposes",
        )
        return

    def run(self, args: argparse.Namespace) -> None:
        """Run subcommand with parsed arguments"""
        # set up logging
        setup_logger(logfile=args.logger, silent=args.silent, debug=args.debug)
        log = get_logger()
        # print information about the run
        from new_majiq._version import version

        log.info(f"new-majiq v{version}")
        log.info(f"Command: {' '.join(sys.argv)}")
        log.info(f"From: {Path(os.getcwd()).resolve()}")
        log.info(
            "\n".join(
                [
                    "Arguments:",
                    "{",
                    *(
                        f" {key} = {value},"
                        for key, value in vars(args).items()
                        if key != "func"
                    ),
                    "}",
                ]
            )
        )
        if getattr(args, "use_dask", False):
            # we are using dask, so set up client here
            client = Client(
                n_workers=1,
                threads_per_worker=args.nthreads,
                dashboard_address=":0",
                memory_limit=args.memory_limit,
                local_directory=args.dask_local_directory,
            )
            # dashboard_address=None should disable server, but it doesn't
            client.cluster.scheduler.http_server.stop()  # stop server
            log.info(client)
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
