import multiprocessing

import gunicorn.app.base
from gunicorn.six import iteritems

from voila.flask_proj import deltapsi


def number_of_workers():
    return (multiprocessing.cpu_count() * 2) + 1


def number_of_threads():
    return multiprocessing.cpu_count() * 2


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


if __name__ == '__main__':
    options = {
        'bind': '%s:%s' % ('127.0.0.1', '8080'),
        'workers': number_of_workers(),
        'threads': number_of_threads(),
        'worker_class': 'gthread'
    }
    StandaloneApplication(deltapsi.app, options).run()
