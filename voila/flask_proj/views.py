from waitress import serve

from voila import constants
from voila.config import Config
from voila.flask_proj import deltapsi, psi


def run_service():
    config = Config()
    analysis_type = config.analysis_type
    if analysis_type == constants.ANALYSIS_PSI:
        serve(psi.app, listen='127.0.0.1:55555', threads=config.nproc)

    elif analysis_type == constants.ANALYSIS_DELTAPSI:
        serve(deltapsi.app, listen='127.0.0.1:55555', threads=config.nproc)
