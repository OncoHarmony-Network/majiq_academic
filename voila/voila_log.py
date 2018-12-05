import logging
from logging import Formatter, StreamHandler, getLogger, Logger
from logging.handlers import RotatingFileHandler
from pathlib import Path

from voila.constants import VOILA_LOG_NAME


def voila_log(filename=None, silent=False, debug=False):
    """
    Logger used throughout voila.  After this has been initialized, then it will retrieve the same logger each time
    this function is called.
    :param debug:
    :param filename: location of log
    :param silent: if true, then logger will not print to command line
    :return: log
    """

    try:
        return Logger.manager.loggerDict[VOILA_LOG_NAME]
    except KeyError:
        pass

    formatter = Formatter("%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s")

    log = getLogger(VOILA_LOG_NAME)

    if filename:
        filename = Path(filename).expanduser().resolve()

        # keep newest 2 gigs of logs in two files
        handler = RotatingFileHandler(filename, maxBytes=1000 * 1000 * 1000, backupCount=2)
        handler.setFormatter(formatter)
        log.addHandler(handler)

    if not silent:
        streamHandler = StreamHandler()
        streamHandler.setFormatter(formatter)
        log.addHandler(streamHandler)

    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    return log
