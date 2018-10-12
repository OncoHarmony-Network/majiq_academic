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


class Config:
    def __init__(self):
        # todo: this should be a singleton class
        self.default = None
        self.analysis_type = None

        self.files = None
        self.voila_files = None
        self.splice_graph_file = None
        self.output = None
        self.nproc = None
        self.threshold = None
        self.non_changing_threshold = None
        self.probability_threshold = None
        self.show_all = None

        self.read()

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

    @staticmethod
    def find_analysis_type(matrix_files):
        analysis_type = None
        for mf in matrix_files:
            with Matrix(mf) as m:
                if analysis_type is None:
                    analysis_type = m.analysis_type
                if analysis_type != m.analysis_type:
                    raise MixedAnalysisTypeVoilaFiles()

        if analysis_type in [constants.ANALYSIS_PSI, constants.ANALYSIS_DELTAPSI]:
            if len(matrix_files) > 1:
                raise FoundMoreThanOneVoilaFile()

        return analysis_type

    @classmethod
    def matrix_files(cls, vs):
        matrix_files = set()

        for v in vs:

            v = Path(v)

            if v.is_file():

                try:
                    with Matrix(v):
                        matrix_files.add(v)
                except OSError:
                    pass

            elif v.is_dir():
                x = cls.matrix_files(v.iterdir())
                matrix_files.update(x)

        return matrix_files

    @classmethod
    def write(cls, args):
        splice_graph_file = cls.splice_graph_file(args.files)
        matrix_files = cls.matrix_files(args.files)
        analysis_type = cls.find_analysis_type(matrix_files)
        config = configparser.ConfigParser()
        files = 'FILES'
        default = 'DEFAULT'

        config.add_section(files)
        config.set(files, 'voila', '\n'.join(str(m) for m in matrix_files))
        config.set(files, 'splice_graph', str(splice_graph_file))

        config.set(default, 'analysis_type', analysis_type)
        config.set(default, 'output', args.output)
        config.set(default, 'nproc', str(args.nproc))
        config.set(default, 'threshold', str(args.threshold))
        config.set(default, 'non_changing_threshold', str(args.non_changing_threshold))
        config.set(default, 'probability_threshold', str(args.probability_threshold))
        config.set(default, 'show_all', str(args.show_all))

        with open(constants.CONFIG_FILE, 'w') as configfile:
            config.write(configfile)

    def read(self):
        config = configparser.ConfigParser()
        config.read(constants.CONFIG_FILE)

        self.files = config['FILES']
        self.voila_files = self.files['voila'].split('\n')
        self.splice_graph_file = self.files['splice_graph']

        self.default = config['DEFAULT']
        self.analysis_type = self.default['analysis_type']
        self.output = self.default['output']
        self.nproc = self.default['nproc']
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
