import configparser
from pathlib import Path

import sqlalchemy
from sqlalchemy.exc import DatabaseError

from voila import constants
from voila.api import Matrix, SpliceGraph


class FoundNoSpliceGraphFile(Exception):
    pass


class FoundMoreThanOneSpliceGraph(Exception):
    pass


class MixedAnalysisTypeVoilaFiles(Exception):
    pass


class FoundMoreThanOneVoilaFile(Exception):
    pass


class Singleton(object):
    def __new__(cls, *args, **kwds):
        it = cls.__dict__.get("__it__")
        if it is not None:
            return it
        cls.__it__ = it = object.__new__(cls)
        it.init(*args, **kwds)
        return it

    def init(self, *args, **kwds):
        pass


class Config(Singleton):
    def __new__(cls, *args, **kwargs):
        config = configparser.ConfigParser()
        config.read(constants.CONFIG_FILE)
        analysis_type = config['DEFAULT']['analysis_type']

        if analysis_type == constants.ANALYSIS_PSI:
            c = super().__new__(PsiConfig)

        elif analysis_type == constants.ANALYSIS_DELTAPSI:
            c = super().__new__(DeltaPsiConfig)

        elif analysis_type == constants.ANALYSIS_HETEROGEN:
            c = super().__new__(HeterogenConfig)
        else:
            raise Exception()

        files = config['FILES']
        c.default = config['DEFAULT']

        c.voila_files = files['voila'].split('\n')
        c.splice_graph_file = files['splice_graph']

        c.analysis_type = c.default['analysis_type']
        c.output = c.default['output']
        c.nproc = int(c.default['nproc'])

        return c

    @property
    def voila_file(self):
        return self.voila_files[0]

    @classmethod
    def splice_graph_file(cls, vs):
        sg_files = set()

        for v in vs:

            v = Path(v)

            if v.is_file():

                try:
                    with SpliceGraph(v):
                        sg_files.add(v)
                except sqlalchemy.exc.DatabaseError:
                    pass

            elif v.is_dir():

                try:
                    v_sg_file = cls.splice_graph_file(v.iterdir())
                    sg_files.add(v_sg_file)
                except FoundNoSpliceGraphFile:
                    pass

        if len(sg_files) == 0:
            raise FoundNoSpliceGraphFile()

        if len(sg_files) > 1:
            raise FoundMoreThanOneSpliceGraph()

        sg_file = sg_files.pop()

        return sg_file.resolve()

    @classmethod
    def find_voila_files(cls, vs):
        voila_files = set()

        for v in vs:

            v = Path(v)

            if v.is_file():

                try:
                    with Matrix(v):
                        voila_files.add(v)
                except OSError:
                    pass

            elif v.is_dir():
                x = cls.find_voila_files(v.iterdir())
                voila_files.update(x)

        return voila_files

    @staticmethod
    def find_analysis_type(voila_files):
        analysis_type = None

        for mf in voila_files:

            with Matrix(mf) as m:

                if analysis_type is None:
                    analysis_type = m.analysis_type

                if analysis_type != m.analysis_type:
                    raise MixedAnalysisTypeVoilaFiles()

        if analysis_type in [constants.ANALYSIS_PSI, constants.ANALYSIS_DELTAPSI]:

            if len(voila_files) > 1:
                raise FoundMoreThanOneVoilaFile()

        return analysis_type

    @classmethod
    def _analysis_type_config(cls, args, config):
        raise NotImplementedError()

    @classmethod
    def write(cls, args):
        voila_files = cls.find_voila_files(args.files)
        analysis_type = cls.find_analysis_type(voila_files)
        splice_graph_file = cls.splice_graph_file(args.files)

        config = configparser.ConfigParser()
        files = 'FILES'
        default = 'DEFAULT'

        config.add_section(files)
        config.set(files, 'voila', '\n'.join(str(m) for m in voila_files))
        config.set(files, 'splice_graph', str(splice_graph_file))

        config.set(default, 'analysis_type', analysis_type)
        config.set(default, 'output', args.output)
        config.set(default, 'nproc', str(args.nproc))

        if analysis_type == constants.ANALYSIS_PSI:
            analysis_type_config = PsiConfig._analysis_type_config

        elif analysis_type == constants.ANALYSIS_DELTAPSI:
            analysis_type_config = DeltaPsiConfig._analysis_type_config

        elif analysis_type == constants.ANALYSIS_HETEROGEN:
            analysis_type_config = HeterogenConfig._analysis_type_config

        else:
            raise Exception()

        analysis_type_config(args, config)

        with open(constants.CONFIG_FILE, 'w') as configfile:
            config.write(configfile)


class PsiConfig(Config):
    @classmethod
    def _analysis_type_config(cls, args, config):
        pass


class DeltaPsiConfig(Config):
    def __init__(self):
        self.threshold = float(self.default['threshold'])
        self.non_changing_threshold = float(self.default['non_changing_threshold'])

        try:
            self.probability_threshold = float(self.default['probability_threshold'])
        except ValueError:
            if self.default['probability_threshold'] == 'None':
                self.probability_threshold = None
            else:
                raise

        self.show_all = self.default['show_all'] == 'True'

    @classmethod
    def _analysis_type_config(cls, args, config):

        default = 'DEFAULT'
        config.set(default, 'threshold', str(args.threshold))
        config.set(default, 'non_changing_threshold', str(args.non_changing_threshold))
        config.set(default, 'probability_threshold', str(args.probability_threshold))
        config.set(default, 'show_all', str(args.show_all))


class HeterogenConfig(Config):
    @classmethod
    def _analysis_type_config(cls, args, config):
        pass
