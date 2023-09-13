import logging
from logging import Formatter, StreamHandler, getLogger, Logger
from logging.handlers import RotatingFileHandler
from pathlib import Path
import os

from rna_voila.constants import VOILA_LOG_NAME

# this is only required to avoid BasicConfig call in the GtfParse library ruining all of our logger configs
# may be removed in the future when presumably this bug would be fixed.
# https://github.com/openvax/gtfparse/issues/40
logging.getLogger().addHandler(logging.NullHandler())

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
    log.setLevel(logging.DEBUG)

    if filename:
        filename = Path(filename).expanduser().resolve()
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        # keep newest 2 gigs of logs in two files
        handler = RotatingFileHandler(filename, maxBytes=1000 * 1000 * 1000, backupCount=2)
        handler.setFormatter(formatter)
        handler.setLevel(logging.DEBUG)
        log.addHandler(handler)

    if not silent:
        streamHandler = StreamHandler()
        streamHandler.setFormatter(formatter)
        if debug:
            streamHandler.setLevel(logging.DEBUG)
        else:
            streamHandler.setLevel(logging.INFO)
        log.addHandler(streamHandler)



    return log
