import configparser
import inspect
import sqlite3
from collections import namedtuple
from pathlib import Path
import sys
from rna_voila import constants

from rna_voila.api import Matrix, SpliceGraph
from rna_voila.api.splice_graph_lr import SpliceGraphLR
from rna_voila.exceptions import FoundNoSpliceGraphFile, FoundMoreThanOneSpliceGraph, \
    MixedAnalysisTypeVoilaFiles, FoundMoreThanOneVoilaFile, AnalysisTypeNotFound
from rna_voila.voila_log import voila_log


_log_keys = ['logger', 'silent']
_sys_keys = ['nproc', 'debug']
_global_keys = ['analysis_type', 'memory_map_hdf5']

_ViewConfig = namedtuple('ViewConfig', _global_keys + _sys_keys + _log_keys + ['voila_file', 'voila_files',
                                        'splice_graph_file',
                                        'force_index', 'port', 'host', 'web_server', 'index_file',
                                        'num_web_workers', 'strict_indexing', 'skip_type_indexing', 'splice_graph_only',
                                        'enable_passcode', 'ignore_inconsistent_group_errors',
                                        'enable_het_comparison_chooser', 'long_read_file',  'disable_reads',
                                        'group_order_override'])
_ViewConfig.__new__.__defaults__ = (None,) * len(_ViewConfig._fields)
_TsvConfig = namedtuple('TsvConfig', _global_keys + _sys_keys + _log_keys + ['file_name', 'voila_files', 'voila_file',
                                      'splice_graph_file',
                                      'non_changing_threshold', 'threshold', 'show_all',
                                      'probability_threshold', 'gene_ids', 'gene_names', 'lsv_ids',
                                      'lsv_types', 'strict_indexing', 'show_read_counts',
                                      'ignore_inconsistent_group_errors', 'non_changing_pvalue_threshold',
                                      'non_changing_within_group_iqr',
                                      'non_changing_between_group_dpsi', 'changing_pvalue_threshold',
                                      'changing_between_group_dpsi'])
_TsvConfig.__new__.__defaults__ = (None,) * len(_TsvConfig._fields)
_ClassifyConfig = namedtuple('ClassifyConfig', _global_keys + _sys_keys + _log_keys + ['directory', 'voila_files',
                                      'voila_file', 'splice_graph_file',
                                      'decomplexify_psi_threshold', 'decomplexify_deltapsi_threshold',
                                      'decomplexify_reads_threshold', 'gene_ids',
                                      'keep_constitutive', 'keep_no_lsvs_modules', 'only_binary',
                                      'untrimmed_exons', 'putative_multi_gene_regions',
                                        'probability_changing_threshold',
                                        'probability_non_changing_threshold', 'show_all',
                                        'non_changing_pvalue_threshold', 'non_changing_within_group_iqr',
                                        'non_changing_between_group_dpsi', 'changing_pvalue_threshold',
                                        'changing_between_group_dpsi', 'changing_between_group_dpsi_secondary',
                                        'keep_no_lsvs_junctions', 'debug_num_genes', 'overwrite', 'output_mpe',
                                        'heatmap_selection', 'enabled_outputs',
                                        'ignore_inconsistent_group_errors', 'disable_metadata',
                                        'show_read_counts', 'cassettes_constitutive_column',
                                        'non_changing_median_reads_threshold', 'permissive_event_non_changing_threshold',
                                        'include_change_cases', 'junc_gene_dist_column'])
_ClassifyConfig.__new__.__defaults__ = (None,) * len(_ClassifyConfig._fields)
_FilterConfig = namedtuple('FilterConfig', _global_keys + _sys_keys + _log_keys + ['directory', 'voila_files',
                                            'voila_file', 'splice_graph_file',
                                            'gene_ids', 'overwrite',
                                            'gene_ids_file', 'lsv_ids', 'lsv_ids_file', 'voila_files_only',
                                            'splice_graph_only',
                                            'changing_threshold', 'non_changing_threshold',
                                            'probability_changing_threshold',
                                            'probability_non_changing_threshold', 'changing', 'non_changing'])
_FilterConfig.__new__.__defaults__ = (None,) * len(_FilterConfig._fields)
_SplitterConfig = namedtuple('SplitterConfig', _global_keys + _sys_keys + _log_keys + ['directory', 'voila_files',
                                        'voila_file', 'splice_graph_file', 'num_divisions', 'copy_only', 'overwrite'])
_SplitterConfig.__new__.__defaults__ = (None,) * len(_SplitterConfig._fields)
_RecombineConfig = namedtuple('RecombineConfig', _global_keys + _sys_keys + _log_keys + ['directories', 'directory'])
_RecombineConfig.__new__.__defaults__ = (None,) * len(_RecombineConfig._fields)
_LongReadsConfig = namedtuple('LongReadsConfig', _global_keys + _sys_keys + _log_keys + ['voila_file', 'lr_gtf_file', 'lr_bed_file', 'lr_tsv_file', 'splice_graph_file', 'output_file', 'gene_id'])
_LongReadsConfig.__new__.__defaults__ = (None,) * len(_LongReadsConfig._fields)

# global config variable to act as the singleton instance of the config.
this_config = None


def find_splice_graph_file(vs):
    """
    Function that located all splice graphs from a list of files and directories.
    :param vs: list of files and directories.
    :return: location of splice graph file
    """

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


def find_voila_files(vs, group_order_override=None):
    """
    Find all voila files in files and directories.
    :param vs: list of files and directories.
    :return: list of voila files
    """

    voila_files = []
    voila_files_to_group_names = {}

    for v in vs:
        v = Path(v)

        if v.is_file() and v.name.endswith('.voila'):

            try:
                with Matrix(v, pre_config=True) as m:
                    voila_files.append(v)
                    if group_order_override:
                        voila_files_to_group_names[v] = m.group_names[0]
            except OSError:
                voila_log().warning('Error opening voila file %s , skipping this file' % str(v))
                pass

        elif v.is_dir():
            x, x2 = find_voila_files(v.iterdir(), group_order_override)
            voila_files = [*voila_files, *x]
            voila_files_to_group_names.update(x2)



    # We rely on the directory of voila files to store the index for het runs, therefore it would be best to
    # have the same directory every time.
    voila_files.sort()

    return voila_files, voila_files_to_group_names

def reorder_voila_files(voila_files, group_order_override, voila_files_to_group_names):
    try:
        return sorted(voila_files, key=lambda x: group_order_override.index(voila_files_to_group_names[x]))
    except ValueError:
        voila_log().critical("Could not match group order override to provided voila files")
        raise


def get_mixed_analysis_type_str(voila_files):
    types = {'psi': 0, 'delta_psi': 0, 'het': 0}
    for mf in voila_files:

        with Matrix(mf, pre_config=True) as m:

            if m.analysis_type == constants.ANALYSIS_PSI:
                types['psi'] += 1

            elif m.analysis_type == constants.ANALYSIS_DELTAPSI:
                types['delta_psi'] += 1

            elif m.analysis_type == constants.ANALYSIS_HETEROGEN:
                types['het'] += 1

    strsout = []
    if types['psi']:
        strsout.append("PSIx%d" % types['psi'])
    if types['delta_psi']:
        strsout.append("dPSIx%d" % types['delta_psi'])
    if types['het']:
        strsout.append("HETx%d" % types['het'])
    return ' '.join(strsout)

def find_analysis_type(voila_files):
    """
    Find the analysis type from the voila files.
    :param voila_files: list of voila files.
    :return: String
    """
    analysis_type = None

    for mf in voila_files:

        with Matrix(mf, pre_config=True) as m:

            if analysis_type is None:
                analysis_type = m.analysis_type

            if analysis_type != m.analysis_type:
                raise MixedAnalysisTypeVoilaFiles()

    if not analysis_type:
        raise AnalysisTypeNotFound()

    if analysis_type in (constants.ANALYSIS_DELTAPSI,) and len(voila_files) > 1:
        raise FoundMoreThanOneVoilaFile()

    return analysis_type


def write(args):
    """
    Write command line argmuments into a ini file.
    :param args: argparse object
    :return: None
    """

    voila_log().info('config file: ' + constants.CONFIG_FILE)

    attrs = inspect.getmembers(args, lambda a: not inspect.isbuiltin(a))
    attrs = (a for a in attrs if not a[0].startswith('_'))
    attrs = dict(attrs)

    if not args.func.__name__ in ('recombine', 'longReadsInputsToLongReadsVoila'):
        sg_file = find_splice_graph_file(args.files)
    else:
        sg_file = None

    if (hasattr(args, 'splice_graph_only') and args.splice_graph_only) or args.func.__name__ in ('recombine', 'longReadsInputsToLongReadsVoila'):

        analysis_type = ''
        voila_files = []

    else:


        group_order_override = getattr(args, "group_order_override", None)
        voila_files, voila_files_to_group_names = find_voila_files(args.files, group_order_override=group_order_override)
        if group_order_override:
            voila_files = reorder_voila_files(voila_files, group_order_override, voila_files_to_group_names)
        if args.func.__name__ in ("Filter", "Classify", 'splitter', 'recombine'):
            analysis_type = get_mixed_analysis_type_str(voila_files)
        else:
            analysis_type = find_analysis_type(voila_files)

    # raise multi-file error if trying to run voila in TSV mode with multiple input files
    # (currently, multiple input is only supported in View mode)
    if analysis_type in (constants.ANALYSIS_PSI, ) and \
            args.func.__name__ not in ['run_service', 'Classify', 'splitter', 'recombine'] and len(voila_files) > 1:

        raise FoundMoreThanOneVoilaFile()

    # attributes that don't need to be in the ini file
    for remove_key in ['files', 'func', 'logger']:
        try:
            del attrs[remove_key]
        except KeyError:
            pass

    config_parser = configparser.ConfigParser()
    files = 'FILES'
    settings = 'SETTINGS'
    filters = 'FILTERS'

    # Get filters from arguments, add them to the appropriate section, and remove them from arguments.
    for lsv_filter in ['lsv_types', 'lsv_ids', 'gene_ids', 'gene_names']:
        if lsv_filter in attrs and attrs[lsv_filter]:
            try:
                config_parser.set(filters, lsv_filter, '\n'.join(attrs[lsv_filter]))
            except configparser.NoSectionError:
                config_parser.add_section(filters)
                config_parser.set(filters, lsv_filter, '\n'.join(attrs[lsv_filter]))

            del attrs[lsv_filter]

    if attrs.get('enable_passcode', False):
        import random, string
        attrs['enable_passcode'] = ''.join(random.choices(string.ascii_letters + string.digits, k=20))
    else:
        attrs['enable_passcode'] = ''

    # Get settings from arguments.
    config_parser.add_section(settings)
    for key, value in attrs.items():
        if isinstance(value, int) or isinstance(value, float) or value:
            config_parser.set(settings, key, str(value))

    if args.func.__name__ == "Classify":
        # check if default = default for some args depending on the type of classify being run.
        # we store a default of none and apply the actual default if the user does not specify another value
        if not config_parser.has_option(settings, 'decomplexify_psi_threshold'):
            if config_parser.getboolean(settings, 'show_all'):
                config_parser.set(settings, 'decomplexify_psi_threshold', '0.0')

        if not config_parser.has_option(settings, 'decomplexify_deltapsi_threshold'):
            if config_parser.getboolean(settings, 'show_all'):
                config_parser.set(settings, 'decomplexify_deltapsi_threshold', '0.0')

    config_parser.set(settings, 'analysis_type', analysis_type)

    # Get files from arguments
    config_parser.add_section(files)
    config_parser.set(files, 'voila', '\n'.join(str(m) for m in voila_files))
    config_parser.set(files, 'splice_graph', str(sg_file))

    # Write ini file.
    with open(constants.CONFIG_FILE, 'w') as configfile:
        config_parser.write(configfile)

    # dump config file to log
    with open(constants.CONFIG_FILE, 'r') as configfile:
        for line in configfile:
            voila_log().debug('CFG| ' + line[:-1])

    if config_parser.has_option(settings, 'long_read_file') and config_parser.get(settings, 'long_read_file'):
        voila_log().info(f"Parsing long reads file: {config_parser.get(settings, 'long_read_file')}")
        SpliceGraphLR(config_parser.get(settings, 'long_read_file'))



class ViewConfig:
    def __new__(cls, *args, **kwargs):
        """
        Before the object is created, we'll parse the ini file, save the named tuple to a global variable, and use it
        as the sington object. This class is specifically for the HTML view.

        :param args: arguments
        :param kwargs: keyword arguments
        :return: named tuple config
        """

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
            for int_key in ['nproc', 'port', 'num_web_workers']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for bool_key in ['force_index', 'silent', 'debug', 'strict_indexing', 'skip_type_indexing',
                             'ignore_inconsistent_group_errors', 'enable_het_comparison_chooser', 'memory_map_hdf5',
                             'disable_reads']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            # singleton data store properties
            if settings.get('memory_map_hdf5', False) and not 'index_file' in settings:
                voila_log().critical('To use hdf5 memory map performance mode, you must specify --index-file as well')
                sys.exit(1)

            this_config = _ViewConfig(**{**files, **settings})

        return this_config


class TsvConfig:
    def __new__(cls, *args, **kwargs):
        """
        Before the object is created, we'll parse the ini file, save the named tuple to a global variable, and use it
        as the sington object. This class is specifically for the TSV output.

        :param args: arguments
        :param kwargs: keyword arguments
        :return: named tuple config
        """

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
            for float_key in ['non_changing_threshold', 'threshold', 'probability_threshold',
                              'non_changing_pvalue_threshold',
                              'non_changing_within_group_iqr', 'non_changing_between_group_dpsi',
                              'changing_pvalue_threshold', 'changing_between_group_dpsi']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['show_all', 'silent', 'debug', 'strict_indexing', 'show_read_counts',
                             'ignore_inconsistent_group_errors', 'memory_map_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _TsvConfig(**{**files, **settings, **filters})

        return this_config

class ClassifyConfig:
    def __new__(cls, *args, **kwargs):
        """

        """

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


            for int_key in ['nproc', 'keep_constitutive', 'decomplexify_reads_threshold', 'debug_num_genes',
                            'non_changing_median_reads_threshold']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in ['decomplexify_psi_threshold', 'decomplexify_deltapsi_threshold',
                              'probability_changing_threshold',
                              'probability_non_changing_threshold', 'non_changing_pvalue_threshold',
                              'non_changing_within_group_iqr', 'non_changing_between_group_dpsi',
                              'changing_pvalue_threshold', 'changing_between_group_dpsi',
                              'changing_between_group_dpsi_secondary', 'permissive_event_non_changing_threshold']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'keep_no_lsvs_modules', 'only_binary', 'untrimmed_exons', 'overwrite',
                             'putative_multi_gene_regions', 'show_all', 'keep_no_lsvs_junctions', 'output_mpe',
                             'ignore_inconsistent_group_errors', 'disable_metadata', 'show_read_counts',
                             'cassettes_constitutive_column', 'include_change_cases', 'junc_gene_dist_column',
                             'memory_map_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            if settings['decomplexify_reads_threshold'] == 0:
                voila_log().warning("--decomplexify-reads-threshold 0 is not recommended and not tested!")

            # implications
            if settings['putative_multi_gene_regions']:
                settings['keep_constitutive'] = True
            if settings['keep_constitutive']:
                settings['keep_no_lsvs_modules'] = True
                settings['keep_no_lsvs_junctions'] = True



            if not settings['putative_multi_gene_regions']:
                settings['enabled_outputs'] = ['summary', 'events', 'junctions', 'heatmap']
                if settings['output_mpe']:
                    settings['enabled_outputs'].append('mpe')
                    settings['keep_constitutive'] = True
                    settings['keep_no_lsvs_modules'] = True
                    settings['keep_no_lsvs_junctions'] = True

            if not settings['show_all']:
                if 'HET' not in settings['analysis_type'] and 'dPSI' not in settings['analysis_type']:
                    voila_log().warning("Only PSI files were provided, so unable to detect changing events. "
                                        "Enabling --show-all automatically. ")
                    settings['show_all'] = True

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _ClassifyConfig(**{**files, **settings, **filters})

        return this_config

class FilterConfig:
    def __new__(cls, *args, **kwargs):
        """
        Before the object is created, we'll parse the ini file, save the named tuple to a global variable, and use it
        as the sington object. This class is specifically for the TSV output.

        :param args: arguments
        :param kwargs: keyword arguments
        :return: named tuple config
        """
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
            for float_key in ['non_changing_threshold', 'changing_threshold', 'probability_changing_threshold',
                              'probability_non_changing_threshold']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'overwrite', 'voila_files_only', 'splice_graph_only', 'changing', 'non_changing',
                             'memory_map_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _FilterConfig(**{**files, **settings, **filters})

        return this_config

class SplitterConfig:
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

            for int_key in ['nproc', 'num_divisions']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in []:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'copy_only', 'overwrite', 'memory_map_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _SplitterConfig(**{**files, **settings, **filters})

        return this_config


class RecombineConfig:
    def __new__(cls, *args, **kwargs):

        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            settings = dict(config_parser['SETTINGS'])

            for int_key in ['nproc']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in []:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'memory_map_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _RecombineConfig(**{**settings, **filters})

        return this_config

class LongReadsConfig:
    def __new__(cls, *args, **kwargs):

        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            settings = dict(config_parser['SETTINGS'])

            for int_key in []:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in []:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in []:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _LongReadsConfig(**{**settings, **filters})

        return this_config