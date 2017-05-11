import errno
import logging
import os
from logging import Formatter, StreamHandler, getLogger, Logger
from logging.handlers import RotatingFileHandler

from voila.constants import VOILA_LOG_NAME


def voila_log(filename=None, silent=False):
    """
    Logger used throughout voila.  After this has been initialized, then it will retrieve the same logger each time
    this function is called.
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
        try:
            log_directory = os.path.dirname(filename)
            if log_directory:
                os.makedirs(log_directory)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        # keep newest 2 gigs of logs in two files
        handler = RotatingFileHandler(filename, maxBytes=1000 * 1000 * 1000, backupCount=2)
        handler.setFormatter(formatter)
        log.addHandler(handler)

    if not silent:
        streamHandler = StreamHandler()
        streamHandler.setFormatter(formatter)
        log.addHandler(streamHandler)

    log.setLevel(logging.DEBUG)

    return log
