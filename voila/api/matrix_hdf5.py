import os
from abc import ABC, abstractmethod
from typing import List

import h5py
import numpy

from voila import constants
from voila.exceptions import GeneIdNotFoundInVoilaFile
from voila.vlsv import collapse_matrix


def lsv_id_to_gene_id(lsv_id):
    return ':'.join(lsv_id.split(':')[:-2])


class MatrixHdf5:
    def __init__(self, filename, mode='r', lock=None):
        """
        Access voila's HDF5 file.

        :param filename: name of voila file
        :param mode: generally r or w
        """
        filename = os.path.expanduser(filename)
        self.dt = h5py.special_dtype(vlen=numpy.unicode)
        self.h = h5py.File(filename, mode, libver='latest')
        self.lock = lock

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.h.close()

    def _add(self, lsv_id: str, key: str, data):
        gene_id = lsv_id_to_gene_id(lsv_id)
        try:
            lsvs_grp = self.h['lsvs']
        except KeyError:
            lsvs_grp = self.h.create_group('lsvs')

        try:
            gene_id_grp = lsvs_grp[gene_id]
        except KeyError:
            gene_id_grp = lsvs_grp.create_group(gene_id)

        try:
            lsv_id_grp = gene_id_grp[lsv_id]
        except KeyError:
            lsv_id_grp = gene_id_grp.create_group(lsv_id)

        lsv_id_grp.create_dataset(key, data=data)

    def add(self, lsv_id: str, key: str, data):
        """
        Add a key/value pair for a LSV ID.

        :param lsv_id: unique LSV identifier
        :param key: key for value
        :param data: data to store
        :return: None
        """
        try:
            self.lock.acquire()
            self._add(lsv_id, key, data)
            self.lock.release()
        except AttributeError:
            self._add(lsv_id, key, data)

    def add_multiple(self, lsv_id, **kwargs):
        """
        Add multiple key/values for a LSV ID.

        :param lsv_id: unique LSV identifier
        :param kwargs: keyword argument storing key/values
        :return: None
        """
        for key, value in kwargs.items():
            self.add(lsv_id, key, value)

    def get(self, lsv_id: str, key: str):
        """
        Retrieve value from file for a LSV ID.

        :param lsv_id: unique LSV identifier
        :param key: key for value
        :return: Value
        """
        gene_id = lsv_id_to_gene_id(lsv_id)
        lsv_grp = self.h['lsvs'][gene_id][lsv_id]
        return lsv_grp[key].value

    def get_many(self, lsv_id: str, keys: List[str]):
        """
        Retrieve many values for a LSV ID.

        :param lsv_id: unique LSV identifier
        :param keys: list of keys for values
        :return: Generator of key, values
        """
        gene_id = lsv_id_to_gene_id(lsv_id)
        lsv_grp = self.h['lsvs'][gene_id][lsv_id]
        for key in keys:
            yield key, lsv_grp[key].value

    @property
    def prior(self):
        return self.h['metadata']['prior'].value

    @prior.setter
    def prior(self, ps):
        self.h.create_dataset('metadata/prior', data=[collapse_matrix(p) for p in ps])

    @property
    def analysis_type(self):
        return self.h['metadata']['analysis_type'].value

    @analysis_type.setter
    def analysis_type(self, a):
        self.h.create_dataset('metadata/analysis_type', data=a)

    @property
    def group_names(self):
        return self.h['metadata']['group_names'].value

    @group_names.setter
    def group_names(self, n):
        self.h.create_dataset('metadata/group_names', data=numpy.array(n, dtype=self.dt))

    @property
    def experiment_names(self):
        return self.h['metadata']['experiment_names'].value

    @experiment_names.setter
    def experiment_names(self, n):
        self.h.create_dataset('metadata/experiment_names', data=numpy.array(n, dtype=self.dt))

    @property
    def stat_names(self):
        return self.h['metadata']['stat_names'].value

    @stat_names.setter
    def stat_names(self, s):
        self.h.create_dataset('metadata/stat_names', data=numpy.array(s, dtype=self.dt))

    @property
    def gene_ids(self):
        yield from self.h['lsvs']

    def lsv_ids(self, gene_ids=None):

        if not gene_ids:
            gene_ids = self.h['lsvs']

        lsvs = self.h['lsvs']

        try:
            for gene_id in gene_ids:
                yield from lsvs[gene_id]
        except KeyError:
            raise GeneIdNotFoundInVoilaFile(gene_id)

    @property
    def file_version(self):
        metadata = self.h['metadata']
        try:
            return metadata['file_version'].value
        except KeyError:
            return -1

    @file_version.setter
    def file_version(self, version):
        self.h.create_dataset('metadata/file_version', data=version)


class MatrixType(ABC):
    @abstractmethod
    def __init__(self, matrix_hdf5, lsv_id, fields):
        """
        The analysis type of data found in matrix file.

        :param matrix_hdf5: matrix HDF5 object
        :param lsv_id: unique LSV identifier
        :param fields: fields allowed to be stored in this file
        """
        self.matrix_hdf5 = matrix_hdf5
        self.lsv_id = lsv_id
        self.fields = fields
        self._lsv_type = None

    def add(self, **kwargs):
        """
        Add key/values to Matrix file.

        :param kwargs: keyword args containing key/values
        :return: None
        """
        if all(k in kwargs for k in self.fields) and len(kwargs) == len(self.fields):
            self.matrix_hdf5.add_multiple(self.lsv_id, **kwargs)
        else:
            wrong_fields = set(kwargs.keys()) - set(self.fields)
            missing_fields = set(self.fields) - set(kwargs.keys()) - wrong_fields
            msg = []
            if wrong_fields:
                msg.append('Wrong field(s): {0}'.format(', '.join(wrong_fields)))
            if missing_fields:
                msg.append('Missing field(s): {0}'.format(', '.join(missing_fields)))

            raise Exception('; '.join(msg))

    def get(self, key: str):
        """
        Retrieve value from Matrix file.

        :param key: Key for value
        :return: Value
        """
        return self.matrix_hdf5.get(self.lsv_id, key)

    def get_many(self, keys: List[str]):
        """
        Retrieve many values using list of keys

        :param keys: list of keys for values
        :return: Generator of key, values
        """
        return self.matrix_hdf5.get_many(self.lsv_id, keys)

    @property
    def lsv_type(self):
        if self._lsv_type is None:
            self._lsv_type = self.get('lsv_type')
        return self._lsv_type

    @property
    def intron_retention(self):
        return 'i' == self.lsv_type[-1]

    @property
    def reference_exon(self):
        for coord in self.lsv_id.split(':')[-1].split('-'):
            if coord == 'nan':
                yield -1
            else:
                yield int(coord)

    @property
    def is_target(self):
        return self.lsv_id.split(':')[-2] == 't'

    @property
    def gene_id(self):
        return lsv_id_to_gene_id(self.lsv_id)

    @property
    def _prime5(self):
        try:
            lsv_type = self.lsv_type[2:]
            if lsv_type[-1] == 'i':
                lsv_type = lsv_type[:-2]

            splice_sites = set(j[0] for j in lsv_type.split('|'))
            return len(splice_sites) > 1
        except IndexError:
            if self.lsv_type == constants.NA_LSV:
                return constants.NA_LSV
            raise

    @property
    def _prime3(self):
        try:
            lsv_type = self.lsv_type[2:]
            if lsv_type[-1] == 'i':
                lsv_type = lsv_type[:-2]

            exons = {}
            for x in lsv_type.split('|'):
                juncs = x[1:].split('.')
                try:
                    exons[juncs[0]].add(juncs[1])
                except KeyError:
                    exons[juncs[0]] = set(juncs[1])

            for value in exons.values():
                if len(value) > 1:
                    return True

            return False
        except IndexError:
            if self.lsv_type == constants.NA_LSV:
                return constants.NA_LSV
            raise

    @property
    def prime5(self):
        if self.target:
            return self._prime5
        else:
            return self._prime3

    @property
    def prime3(self):
        if self.target:
            return self._prime3
        else:
            return self._prime5

    @property
    def target(self):
        return self.lsv_type[0] == 't'

    @property
    def exon_skipping(self):
        try:
            return self.exon_count > 2
        except TypeError:
            if self.exon_count == constants.NA_LSV:
                return constants.NA_LSV
            raise

    @property
    def exon_count(self):
        try:
            lsv_type = self.lsv_type[2:]
            if lsv_type[-1] == 'i':
                lsv_type = lsv_type[:-2]

            return len(set(x[1:3] for x in lsv_type.split('|'))) + 1
        except IndexError:
            if self.lsv_type == constants.NA_LSV:
                return constants.NA_LSV
            raise


class DeltaPsi(MatrixHdf5):
    class _DeltaPsi(MatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            """
            Store Delta PSI data.

            :param matrix_hdf5: matrix HDF5 object
            :param lsv_id: unique LSV identifier
            """
            fields = ('bins', 'group_bins', 'group_means', 'lsv_type', 'junctions')
            super().__init__(matrix_hdf5, lsv_id, fields)

        def add(self, **kwargs):
            """
            Add keyword arguments specific for Delta PSI data.

            :param kwargs: keyword arguments
            :return: None
            """
            bins_list = kwargs.get('bins')
            bins = [collapse_matrix(bins) for bins in bins_list]
            kwargs['bins'] = bins
            super().add(**kwargs)

    def delta_psi(self, lsv_id):
        """
        Accessor function for Delta PSI file.

        :param lsv_id:
        :return:
        """
        return self._DeltaPsi(self, lsv_id)


class Psi(MatrixHdf5):
    class _Psi(MatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            """
            Store PSI data.

            :param matrix_hdf5: matrix HDF5 object
            :param lsv_id: unique LSV identifier
            """
            fields = ('bins', 'means', 'lsv_type', 'junctions')
            super().__init__(matrix_hdf5, lsv_id, fields)

    def psi(self, lsv_id):
        """
        Accessor function for PSI file.

        :param lsv_id:
        :return:
        """
        return self._Psi(self, lsv_id)


class Heterogen(MatrixHdf5):
    class _Heterogen(MatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            """
            Store PSI data.

            :param matrix_hdf5: matrix HDF5 object
            :param lsv_id: unique LSV identifier
            """
            fields = ('lsv_type', 'junction_stats', 'groups')
            super().__init__(matrix_hdf5, lsv_id, fields)

    def heterogen(self, lsv_id):
        """
        Accessor function for PSI file.

        :param lsv_id:
        :return:
        """
        return self._Heterogen(self, lsv_id)
