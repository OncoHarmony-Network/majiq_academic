import configparser
import inspect
import sqlite3
from collections import namedtuple
from pathlib import Path

from voila import constants
from voila.api import Matrix, SpliceGraph
from voila.exceptions import FoundNoSpliceGraphFile, FoundMoreThanOneSpliceGraph, \
    MixedAnalysisTypeVoilaFiles, FoundMoreThanOneVoilaFile, AnalysisTypeNotFound
from voila.utils.voila_log import voila_log

_ViewConfig = namedtuple('ViewConfig', ['voila_file', 'voila_files', 'splice_graph_file', 'analysis_type', 'nproc',
                                        'force_index', 'debug', 'silent', 'port'])
_ViewConfig.__new__.__defaults__ = (None,) * len(_ViewConfig._fields)
_TsvConfig = namedtuple('TsvConfig', ['file_name', 'voila_files', 'voila_file', 'splice_graph_file',
                                      'non_changing_threshold', 'nproc', 'threshold', 'analysis_type', 'show_all',
                                      'debug', 'probability_threshold', 'silent', 'gene_ids', 'gene_names', 'lsv_ids',
                                      'lsv_types'])
_TsvConfig.__new__.__defaults__ = (None,) * len(_TsvConfig._fields)

this_config = None


def find_splice_graph_file(vs):
    sg_files = set()

    for v in vs:

        v = Path(v)

        if v.is_file():

            try:
                with SpliceGraph(v):
                    sg_files.add(v)
            except sqlite3.DatabaseError:
                pass

        elif v.is_dir():

            try:
                v_sg_file = find_splice_graph_file(v.iterdir())
                sg_files.add(v_sg_file)
            except FoundNoSpliceGraphFile:
                pass

    if len(sg_files) == 0:
        raise FoundNoSpliceGraphFile()

    if len(sg_files) > 1:
        raise FoundMoreThanOneSpliceGraph()

    sg_file = sg_files.pop()

    return sg_file.resolve()


def find_voila_files(vs):
    voila_files = []

    for v in vs:

        v = Path(v)

        if v.is_file() and v.name.endswith('.voila'):

            try:
                with Matrix(v):
                    voila_files.append(v)
            except OSError:
                pass

        elif v.is_dir():
            x = find_voila_files(v.iterdir())
            voila_files = [*voila_files, *x]

    # We rely on the directory of voila files to store the index for het runs, therefore it would be best to
    # have the same directory every time.
    voila_files.sort()

    return voila_files


def find_analysis_type(voila_files):
    analysis_type = None

    for mf in voila_files:

        with Matrix(mf) as m:

            if analysis_type is None:
                analysis_type = m.analysis_type

            if analysis_type != m.analysis_type:
                raise MixedAnalysisTypeVoilaFiles()

    if not analysis_type:
        raise AnalysisTypeNotFound()

    if analysis_type in [constants.ANALYSIS_PSI, constants.ANALYSIS_DELTAPSI]:

        if len(voila_files) > 1:
            raise FoundMoreThanOneVoilaFile()

    return analysis_type


def write(args):
    voila_log().info('config file: ' + constants.CONFIG_FILE)
    attrs = inspect.getmembers(args, lambda a: not inspect.isbuiltin(a))
    attrs = (a for a in attrs if not a[0].startswith('_'))
    attrs = dict(attrs)

    sg_file = find_splice_graph_file(args.files)

    if hasattr(args, 'splice_graph_only') and args.splice_graph_only:
        analysis_type = ''
        voila_files = []

    else:
        voila_files = find_voila_files(args.files)
        analysis_type = find_analysis_type(voila_files)

    for remove_key in ['files', 'func', 'logger', 'splice_graph_only']:
        try:
            del attrs[remove_key]
        except KeyError:
            pass

    config_parser = configparser.ConfigParser()
    files = 'FILES'
    settings = 'SETTINGS'
    filters = 'FILTERS'

    for filter in ['lsv_types', 'lsv_ids', 'gene_ids', 'gene_names']:
        if filter in attrs and attrs[filter]:
            try:
                config_parser.set(filters, filter, '\n'.join(attrs[filter]))
            except configparser.NoSectionError:
                config_parser.add_section(filters)
                config_parser.set(filters, filter, '\n'.join(attrs[filter]))

            del attrs[filter]

    config_parser.add_section(settings)
    for key, value in attrs.items():
        if isinstance(value, int) or value:
            config_parser.set(settings, key, str(value))
    config_parser.set(settings, 'analysis_type', analysis_type)

    config_parser.add_section(files)
    config_parser.set(files, 'voila', '\n'.join(str(m) for m in voila_files))
    config_parser.set(files, 'splice_graph', str(sg_file))

    with open(constants.CONFIG_FILE, 'w') as configfile:
        config_parser.write(configfile)


class ViewConfig:
    def __new__(cls, *args, **kwargs):
        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            files = {
                'voila_files': config_parser['FILES']['voila'].split('\n'),
                'voila_file': config_parser['FILES']['voila'].split('\n')[0],
                'splice_graph_file': config_parser['FILES']['splice_graph']
            }

            settings = dict(config_parser['SETTINGS'])
            for int_key in ['nproc', 'port']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for bool_key in ['force_index']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            this_config = _ViewConfig(**{**files, **settings})

        return this_config


class TsvConfig:
    def __new__(cls, *args, **kwargs):
        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            files = {
                'voila_files': config_parser['FILES']['voila'].split('\n'),
                'voila_file': config_parser['FILES']['voila'].split('\n')[0],
                'splice_graph_file': config_parser['FILES']['splice_graph']
            }

            settings = dict(config_parser['SETTINGS'])
            for int_key in ['nproc']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in ['non_changing_threshold', 'threshold', 'probability_threshold']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _TsvConfig(**{**files, **settings, **filters})

        return this_config