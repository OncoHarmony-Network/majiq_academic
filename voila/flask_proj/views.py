import gunicorn.app.base
from gunicorn.six import iteritems
from waitress import serve

from voila import constants
from voila.config import Config
from voila.flask_proj import deltapsi, heterogen, psi


def run_service():
    config = Config()
    analysis_type = config.analysis_type

    # kwargs = {
    #     'listen': '127.0.0.1:' + str(config.port),
    #     'threads': config.nproc
    # }
    #
    # if analysis_type == constants.ANALYSIS_PSI:
    #     serve(psi.app, **kwargs)
    #
    # elif analysis_type == constants.ANALYSIS_DELTAPSI:
    #     serve(deltapsi.app, **kwargs)
    #
    # elif analysis_type == constants.ANALYSIS_HETEROGEN:
    #     serve(heterogen.app, **kwargs)

    options = {
        'bind': '%s:%s' % ('127.0.0.1', config.port),
        'workers': number_of_workers(),
        'threads': number_of_threads(),
        'worker_class': 'gthread'
    }

    if analysis_type == constants.ANALYSIS_PSI:
        StandaloneApplication(psi.app, options).run()

    elif analysis_type == constants.ANALYSIS_DELTAPSI:
        StandaloneApplication(deltapsi.app, options).run()

    elif analysis_type == constants.ANALYSIS_HETEROGEN:
        StandaloneApplication(heterogen.app, options).run()


def number_of_workers():
    return (Config().nproc * 2) + 1


def number_of_threads():
    return Config().nproc * 2


class StandaloneApplication(gunicorn.app.base.BaseApplication):

    def init(self, parser, opts, args):
        raise NotImplementedError()

    def __init__(self, application, options=None):
        self.options = options or {}
        self.application = application
        super(StandaloneApplication, self).__init__()

    def load_config(self):
        config = dict([(key, value) for key, value in iteritems(self.options)
                       if key in self.cfg.settings and value is not None])
        for key, value in iteritems(config):
            self.cfg.set(key.lower(), value)

    def load(self):
        return self.application
