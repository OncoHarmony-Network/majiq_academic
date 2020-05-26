"""
logger.py

Functions for logging progress in MAJIQ
"""

import logging
import os
import resource
import sys


def monitor(msg):
    print(
        "MONITOR", msg, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000, "MB"
    )
    sys.stdout.flush()


def create_if_not_exists(my_dir, logger=False):
    """Create a directory path if it does not exist"""
    try:
        if logger:
            logger.debug("\nCreating directory %s..." % my_dir)
        os.makedirs(my_dir)
    except OSError:
        if logger:
            logger.debug("\nDirectory %s already exists..." % my_dir)


def get_logger(logger_name, silent=False, debug=False):
    """
    Returns a logger instance. silent=False will silence the logger, debug will
    give more information intended for debugging purposes.
    """
    logging_format = "%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s"
    logger = logging.getLogger(logger_name)
    formatter = logging.Formatter(logging_format)

    file_handler = logging.FileHandler(logger_name, mode="a")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    if debug:
        logger.setLevel(logging.DEBUG)
    elif silent:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def close_logger(logger):
    handlers = logger.handlers[:]
    for handler in handlers:
        handler.close()
        logger.removeHandler(handler)
    logging.shutdown()
