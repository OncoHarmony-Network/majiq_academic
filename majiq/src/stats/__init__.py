from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
all_stats = [basename(f)[:-3] for f in modules if isfile(f) if basename(f)[:-3] != '__init__']
operator = {}


class StatsFactory:

    factories = {}

    @staticmethod
    def add_factory(id_stat, stats_factory):
        StatsFactory.factories.put[id] = stats_factory

    # A Template Method:
    @staticmethod
    def create_stats(id_stat):
        if not StatsFactory.factories.has_key(id_stat):
            StatsFactory.factories[id_stat] = eval(id_stat + '.Factory()')
        return StatsFactory.factories[id_stat].create()


class Stats(object):
    pass

