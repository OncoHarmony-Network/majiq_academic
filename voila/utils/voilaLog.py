import errno
import logging
import os
from logging import Formatter, StreamHandler, getLogger, Logger
from logging.handlers import RotatingFileHandler

from voila.constants import VOILA_LOG_NAME


def voilaLog(filename=None, level=logging.DEBUG):
    try:
        # noinspection PyUnresolvedReferences
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
        handler = RotatingFileHandler(filename, maxBytes=1000 * 1000 * 1000, backupCount=5)
        handler.setFormatter(formatter)
        log.addHandler(handler)

    streamHandler = StreamHandler()

    streamHandler.setFormatter(formatter)
    log.addHandler(streamHandler)
    log.setLevel(level)

    return log
