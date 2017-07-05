import os
from math import ceil

import h5py
import numpy

from voila import constants
from voila.constants import EXPERIMENT_NAMES
from voila.hdf5 import BinsDataSet, HDF5
from voila.splice_graphics import GeneGraphic
from voila.utils.exceptions import GeneIdNotFoundInVoilaFile, GeneIdNotFoundInSpliceGraphFile, InValidAnalysisType
from voila.utils.voila_log import voila_log
from voila.vlsv import VoilaLsv, get_expected_dpsi


class Voila(object):
    VERSION = '/voila_file_version'
    LSVS = 'lsvs'
    ANALYSIS_TYPE = '/analysis_type'

    def __init__(self, voila_file_name, mode):
        """
        Parse or edit the voila (quantifier output) file.
        :param voila_file_name: location of voila file
        :param mode: file mode passed to h5py
        """
        super(Voila, self).__init__()
        self.mode = mode
        self.file_name = voila_file_name
        self.hdf5 = None
        self.file_version = None
        self.hdf5_lsvs = None

    def __enter__(self):
        self.hdf5 = h5py.File(self.file_name, self.mode)

        if self.VERSION not in self.hdf5:
            if self.mode == constants.FILE_MODE.write:
                self.hdf5[self.VERSION] = constants.VOILA_FILE_VERSION

        try:
            self.file_version = self.hdf5[self.VERSION].value
        except KeyError:
            pass

        return self

    def __exit__(self, type, value, traceback):
        """
        Closes file when program exits with block.
        :param type: 
        :param value: 
        :param traceback: 
        :return: 
        """
        self.close()

    def _metainfo(self):
        metainfo = '/metainfo'
        try:
            return self.hdf5[metainfo]
        except KeyError:
            return self.hdf5.create_group(metainfo)

    def add_stat_names(self, stat_names):
        HDF5.create(self._metainfo().attrs, 'stat_names', stat_names)

    def close(self):
        try:
            self.hdf5.close()
        except ValueError:
            pass

    def add_lsv(self, voilaLsv):
        """
        Add VoilaLsv to Voila file.
        :param voilaLsv: VoilaLsv object
        :return: None
        """
        voilaLsv.to_hdf5(self.hdf5)

    def add_genome(self, genome):
        self._metainfo().attrs['genome'] = genome

    def add_metainfo(self, genome, group1, experiments1, group2=None, experiments2=None):
        """
        Add metainfo to Voila file.
        :param genome: genome where the genes are found
        :param group1: first group name (used in psi and deltapsi)
        :param experiments1: first list of experiment names (used in psi and deltapsi)
        :param group2: second group (only deltapsi)
        :param experiments2: second list of experiment names (only deltapsi)
        :return: None
        """

        h = self.hdf5.create_group('/metainfo')

        h.attrs['genome'] = genome

        h.attrs['group1'] = group1
        experiments1_grp = h.create_group('experiments1')
        for index, item in enumerate(experiments1):
            experiments1_grp.attrs[str(index)] = item

        if group2 and experiments2:
            h.attrs['group2'] = group2
            experiments2_grp = h.create_group('experiments2')
            for index, item in enumerate(experiments2):
                experiments2_grp.attrs[str(index)] = item

    def get_metainfo(self):
        """
        Get metainfo from voila file.
        :return: dict
        """
        voila_log().info('Getting Voila Metainfo from {0} ...'.format(self.file_name))
        return {key: self._metainfo().attrs[key] for key in self._metainfo().attrs.keys()}

    def get_voila_lsv(self, gene_id, lsv_id):
        """
        Get LSV by LSV id.
        :param lsv_id:
        :return: VoilaLsv
        """
        return VoilaLsv.easy_from_hdf5(self.hdf5[self.LSVS][gene_id][lsv_id])

    def get_lsv_count(self, args):
        return sum(1 for _ in self.get_lsvs(args))

    def get_voila_lsvs(self, args=None, gene_id=None):
        return (self.get_voila_lsv(gene_id, lsv_id) for gene_id, lsv_id in self.get_lsvs(args, gene_id))

    def get_gene_ids(self, args=None):
        if args:
            if args.gene_ids:
                return args.gene_ids
            elif args.lsv_ids:
                return (lsv_id.split(':')[0] for lsv_id in args.lsv_ids)

        return self.hdf5[self.LSVS].keys()

    def get_lsvs(self, args=None, gene_id=None):
        """
        Get list of LSVs from voila file.
        :return: list
        """
        lsv_ids, threshold = None, None

        if args:
            lsv_ids = args.lsv_ids

            if hasattr(args, 'show_all') and not args.show_all:
                threshold = args.threshold

        if gene_id:
            gene_ids = (gene_id,)
        else:
            gene_ids = self.get_gene_ids(args)

        for gene_id in gene_ids:
            for lsv_id in self.get_lsv_ids(gene_id):
                # Search for LSV
                if not lsv_ids or lsv_id in lsv_ids:
                    if not threshold or VoilaLsv.is_lsv_changing(self.get_lsv_means(gene_id, lsv_id), threshold):
                        yield gene_id, lsv_id

    def get_lsv_means(self, gene_id, lsv_id):
        trunc_bins = BinsDataSet(self.hdf5[self.LSVS][gene_id][lsv_id]).decode_list()
        bins = VoilaLsv._extend_bins(trunc_bins)
        return tuple(get_expected_dpsi(b) for b in bins)

    def get_lsv(self, gene_id, lsv_id):
        return self.hdf5[self.LSVS][gene_id][lsv_id]

    def get_lsv_type(self, gene_id, lsv_id):
        return self.get_lsv(gene_id, lsv_id).attrs['lsv_type']

    def get_gene_name(self, gene_id, lsv_id):
        return self.get_lsv(gene_id, lsv_id).attrs['name']

    def get_lsv_ids(self, gene_id):
        lsvs = self.hdf5[self.LSVS]
        try:
            return lsvs[gene_id].keys()
        except KeyError:
            raise GeneIdNotFoundInVoilaFile(gene_id)

    def add_experiments(self, group_name, experiment_names):
        m = self._metainfo()
        try:
            HDF5.create(m.attrs, 'group_names', numpy.concatenate((m.attrs['group_names'], [group_name])))
            HDF5.create(m.attrs, 'experiment_names', numpy.concatenate((m.attrs['experiment_names'], [experiment_names])))
        except KeyError:
            HDF5.create(m.attrs, 'group_names', numpy.array([group_name]))
            HDF5.create(m.attrs, 'experiment_names', numpy.array([experiment_names]))

    def set_analysis_type(self, analysis_type):
        if analysis_type in (constants.ANALYSIS_DELTAPSI, constants.ANALYSIS_PSI, constants.ANALYSIS_HETEROGEN):
            self.hdf5[self.ANALYSIS_TYPE] = analysis_type
        else:
            raise InValidAnalysisType()

    def check_version(self):
        if self.file_version != constants.VOILA_FILE_VERSION:
            voila_log().warning('Voila file version isn\'t current.  This will probably cause significant '
                                'issues with the voila output.  It would be best to run quantifier again with the '
                                'current version of MAJIQ.')


class SpliceGraphs(object):
    GENES = '/genes'
    ROOT = '/'
    VERSION = '/splice_graph_file_version'

    def __init__(self, splice_graph_file_name, mode):
        """
        Class for creating and accessing the splice graph file.
        :param splice_graph_file_name: path to splice graph file
        :param mode: mode to pass to hdf5
        """
        super(SpliceGraphs, self).__init__()
        self.file_name = splice_graph_file_name
        self.mode = mode
        self.hdf5 = None
        self.limit = None
        self.gene_ids = None
        self.file_version = None

    def __enter__(self):
        """
        Open hdf5 in with block.
        :return: self
        """
        self.hdf5 = h5py.File(self.file_name, self.mode)

        if self.VERSION not in self.hdf5:
            if self.mode == constants.FILE_MODE.write:
                self.hdf5[self.VERSION] = constants.SPLICE_GRAPH_FILE_VERSION

        try:
            self.file_version = self.hdf5[self.VERSION].value
        except KeyError:
            pass

        return self

    def __exit__(self, type, value, traceback):
        """
        Close when with block exits.
        :param type: unused
        :param value: unused
        :param traceback: unused
        :return: None
        """
        self.close()

    def close(self):
        """
        Close hdf5 file.
        :return: None
        """
        try:
            self.hdf5.close()
        except Exception:
            pass

    def erase_splice_graph_file(self):
        """
        Remove splice graph file and reopen it.
        :return:
        """
        os.remove(self.file_name)
        self.__enter__()

    def add_gene(self, gene):
        """
        Add gene object to splice graph file.
        :param gene: GeneGraphic object
        :return: None
        """
        gene.to_hdf5(self.hdf5)

    def get_page_count(self, args):
        gene_count = 0
        log = voila_log()

        log.debug('Start page count')

        if hasattr(args, 'voila_file'):
            with Voila(args.voila_file, 'r') as v:
                for gene_id in self.get_gene_ids(args):
                    try:
                        if any(v.get_lsvs(args, gene_id)):
                            gene_count += 1
                    except GeneIdNotFoundInVoilaFile:
                        pass

        else:
            log.debug('Gene limit is set to {0}'.format(args.limit))
            for _ in self.get_gene_ids(args):
                gene_count += 1
                if gene_count == args.limit:
                    break

        log.debug('End page count')

        return int(ceil(gene_count / float(constants.MAX_GENES)))

    def get_gene_name(self, gene_id):
        return self.hdf5[self.GENES][gene_id].attrs['name']

    def get_gene_ids(self, args=None):
        if args and args.gene_ids:
            return args.gene_ids

        if args and hasattr(args, 'lsv_ids') and args.lsv_ids:
            return (lsv_id.split(':')[0] for lsv_id in args.lsv_ids)

        return self.hdf5[self.GENES].keys()

    def get_genes(self):
        for gene_id in self.get_gene_ids():
            yield self.get_gene(gene_id)

    def get_paginated_genes(self, args):
        log = voila_log()
        log.debug('Getting paginated genes')

        gene_list = []
        gene_count = 0

        for gene_id in self.get_gene_ids(args):
            log.debug('Found {0}'.format(gene_id))
            gene_list.append(self.get_gene(gene_id))
            gene_count += 1

            if gene_count == args.limit:
                break

            if len(gene_list) == constants.MAX_GENES:
                yield gene_list
                gene_list = []

        if gene_list:
            yield gene_list

    def get_paginated_genes_with_lsvs(self, args):
        log = voila_log()
        log.debug('Getting paginated genes with LSVs')

        gene_list = []
        lsv_dict = {}

        with Voila(args.voila_file, 'r') as v:
            for gene_id in self.get_gene_ids(args):
                try:
                    lsvs = tuple(v.get_lsvs(args, gene_id=gene_id))
                except GeneIdNotFoundInVoilaFile:
                    lsvs = None

                if lsvs:
                    gene = self.get_gene(gene_id)
                    lsv_dict[gene_id] = tuple(v.get_voila_lsv(gene_id, lsv_id) for gene_id, lsv_id in lsvs)
                    gene_list.append(gene)

                if len(gene_list) == constants.MAX_GENES:
                    yield lsv_dict, gene_list
                    gene_list = []
                    lsv_dict = {}

            if gene_list:
                yield lsv_dict, gene_list

    def get_gene(self, gene_id):
        """
        Get gene by its gene id.
        :param gene_id: unique gene id
        :return: GeneGraphics
        """
        genes = self.hdf5[self.GENES]
        try:
            gene = genes[gene_id]
        except KeyError:
            raise GeneIdNotFoundInSpliceGraphFile(gene_id)

        return GeneGraphic.easy_from_hdf5(gene)

    def add_experiment_names(self, experiment_names):
        """
        Add experiment names to splice graph.
        :param experiment_names: list of experiment names
        :return: None
        """
        # self.hdf5[self.ROOT].attrs[EXPERIMENT_NAMES] = list(experiment_names)
        HDF5.create(self.hdf5[self.ROOT].attrs, EXPERIMENT_NAMES, experiment_names)

    def get_experiments(self):
        """
        Get list of experiment names from splice graph.
        :return: list
        """
        return self.hdf5[self.ROOT].attrs[EXPERIMENT_NAMES]

    def check_version(self):
        if self.file_version != constants.SPLICE_GRAPH_FILE_VERSION:
            voila_log().warning('Splice graph file version isn\'t current.  This will probably cause significant '
                                'issues with the voila output.  It would be best to run build again with the current '
                                'version of MAJIQ.')
