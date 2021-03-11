"""
logger.py

Wrapper around Python logging for new-majiq

Author: Joseph K Aicher
"""

import logging
from pathlib import Path
from typing import Optional, Union


def get_logger(allow_setup: bool = True) -> logging.Logger:
    """ Get new-majiq logger

    allow_setup: bool
        Set up logger if it hasn't already been when True
        (set to False when running setup)
    """
    # get logger, which we will add handles to
    logger = logging.getLogger(name="new-majiq")
    if allow_setup and not logger.handlers:
        # if no handlers, then it hasn't been setup yet -- make it silent
        setup_logger(silent=True)
    return logger


def _default_formatter() -> logging.Formatter:
    """ How we want log messsages to look
    """
    # set up formatting
    log_format = "%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s"
    formatter = logging.Formatter(log_format)
    return formatter


def setup_logger(
    logfile: Optional[Union[str, Path]] = None,
    silent: bool = False,
    debug: bool = False,
) -> None:
    """ Setup logger to print output to logfile (if None, stderr, as default)
    """
    # set up handler
    if logfile is None:
        # use stderr
        handler = logging.StreamHandler()
    else:
        handler = logging.FileHandler(logfile, mode="a")
    handler.setFormatter(_default_formatter())
    # get logger, remove existing handlers, add this handler to it
    logger = get_logger(allow_setup=False)
    for old_handler in logger.handlers:
        logger.removeHandler(old_handler)
    logger.addHandler(handler)
    # set up level for logging
    if debug:
        logger.setLevel(logging.DEBUG)
    elif silent:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
    # done with setup
    return None
