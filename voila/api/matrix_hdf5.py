import os
from abc import ABC, abstractmethod
from typing import List

import h5py
import numpy as np

from voila import constants
from voila.exceptions import LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile
from voila.vlsv import collapse_matrix


def lsv_id_to_gene_id(lsv_id):
    return ':'.join(lsv_id.split(':')[:-2])


class MatrixHdf5:
    LSVS = 'lsvs'

    def __init__(self, filename, mode='r'):
        """
        Access voila's HDF5 file.

        :param filename: name of voila file
        :param mode: generally r or w
        """
        filename = os.path.expanduser(filename)
        self.dt = h5py.special_dtype(vlen=np.unicode)
        self.h = h5py.File(filename, mode, libver='latest')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.h.close()

    def add_dataset(self, *args, **kwargs):
        grp = self.h
        for k in args[:-1]:
            try:
                grp = grp[k]
            except KeyError:
                grp = grp.create_group(k)
        grp.create_dataset(args[-1], **kwargs)

    def add(self, lsv_id: str, key: str, data):
        """
        Add a key/value pair for a LSV ID.

        :param lsv_id: unique LSV identifier
        :param key: key for value
        :param data: data to store
        :return: None
        """
        gene_id = lsv_id_to_gene_id(lsv_id)
        self.add_dataset(self.LSVS, gene_id, lsv_id, key, data=data)

    def add_multiple(self, lsv_id, **kwargs):
        """
        Add multiple key/values for a LSV ID.

        :param lsv_id: unique LSV identifier
        :param kwargs: keyword argument storing key/values
        :return: None
        """
        for key, value in kwargs.items():
            try:
                self.add(lsv_id, key, value)
            except TypeError:
                print(key, value)
                raise

    def get(self, lsv_id: str, key: str):
        """
        Retrieve value from file for a LSV ID.

        :param lsv_id: unique LSV identifier
        :param key: key for value
        :return: Value
        """
        gene_id = lsv_id_to_gene_id(lsv_id)

        try:
            gene_grp = self.h['lsvs'][gene_id]
        except KeyError:
            raise GeneIdNotFoundInVoilaFile(self.h.filename, gene_id)

        try:
            lsv_grp = gene_grp[lsv_id]
        except KeyError:
            raise LsvIdNotFoundInVoilaFile(self.h.filename, lsv_id)

        try:
            return lsv_grp[key].value
        except KeyError:
            print(dict(lsv_grp))
            raise

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
        """
        Assumed that this is received in log space.
        :param ps:
        :return:
        """
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
        self.h.create_dataset('metadata/group_names', data=np.array(n, dtype=self.dt))

    @property
    def experiment_names(self):
        return self.h['metadata']['experiment_names'].value

    @experiment_names.setter
    def experiment_names(self, ns):
        ns = [[e.decode('utf-8') for e in es] for es in ns]
        arr = np.empty((2, max(len(n) for n in ns)), dtype=self.dt)
        arr.fill('')
        for i, n in enumerate(ns):
            arr[i][0:len(n)] = n
        self.h.create_dataset('metadata/experiment_names', data=arr)

    @property
    def stat_names(self):
        return self.h['metadata']['stat_names'].value

    @stat_names.setter
    def stat_names(self, s):
        self.h.create_dataset('metadata/stat_names', data=np.array(s, dtype=self.dt))

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
        except (KeyError, ValueError):
            return ()

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
        self.matrix_hdf5.add_multiple(self.lsv_id, **kwargs)

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
        coords = self.lsv_id.split(':')[-1].split('-')
        if len(coords) == 2:
            for coord in coords:
                if coord == 'nan':
                    yield -1
                else:
                    yield int(coord)
        else:
            for coord in coords:
                if coord:
                    if coord == '1':
                        yield -1
                    else:
                        yield int(coord)

    @property
    def target(self):
        return self.lsv_id.split(':')[-2] == 't'

    @property
    def gene_id(self):
        return lsv_id_to_gene_id(self.lsv_id)

    @property
    def a5ss(self):
        return 'A5SS' in [self.reference_exon_ss(), self.other_exons_ss()]

    @property
    def a3ss(self):
        return 'A3SS' in [self.reference_exon_ss(), self.other_exons_ss()]

    def reference_exon_ss(self):
        try:
            ss = filter(lambda x: x != 'i', self.lsv_type.split('|')[1:])
            ss = map(lambda x: x.split('.')[0].split('e')[0], ss)
            if len(set(ss)) > 1:
                if self.lsv_type[0] == 's':
                    return 'A5SS'
                else:
                    return 'A3SS'
        except IndexError:
            if self.lsv_type == constants.NA_LSV:
                return constants.NA_LSV
            raise

    def other_exons_ss(self):
        try:
            ss = filter(lambda x: x != 'i', self.lsv_type.split('|')[1:])
            exons = {}
            for x in ss:
                exon = x.split('.')[0].split('e')[1]
                ss = x.split('.')[1]
                try:
                    exons[exon].add(ss)
                except KeyError:
                    exons[exon] = {ss}

            if any(len(values) > 1 for values in exons.values()):
                if self.lsv_type[0] == 's':
                    return 'A3SS'
                else:
                    return 'A5SS'
        except IndexError:
            if self.lsv_type == constants.NA_LSV:
                return constants.NA_LSV
            raise

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
            exons = filter(lambda x: x != 'i', self.lsv_type.split('|')[1:])
            exons = map(lambda x: x.split('.')[0].split('e')[1], exons)
            return len(set(exons)) + 1
        except IndexError:
            if self.lsv_type == constants.NA_LSV:
                return constants.NA_LSV
            raise

    @property
    def binary(self):
        return len(self.lsv_type.split('|')[1:]) == 2


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

        # def add(self, **kwargs):
        #     """
        #     Add keyword arguments specific for Delta PSI data.
        #
        #     :param kwargs: keyword arguments
        #     :return: None
        #     """
        #     bins_list = kwargs.get('bins')
        #     # bins = [collapse_matrix(bins) for bins in bins_list]
        #     # kwargs['bins'] = bins
        #     kwargs['bins'] = bins_list
        #     super().add(**kwargs)

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
            fields = ('lsv_type', 'junction_stats', 'mu_psi', 'mean_psi', 'junctions')
            super().__init__(matrix_hdf5, lsv_id, fields)

        # def add(self, **kwargs):
        #     """
        #     mu_psi: numpy array of two lists. one list for each condition.
        #     junction_stats: 2d numpy matrix. Columns are stat flask_proj values and a row for each junction.
        #     mean_psi: 2d matrix
        #
        #     :param kwargs:
        #     :return:
        #     """
        #
        #     mu_psi = kwargs.get('mu_psi', None)
        #     if mu_psi is not None:
        #         arr = np.empty((2, max(x.shape[0] for x in mu_psi), mu_psi[0].shape[1]))
        #         arr.fill(-1)
        #         for i, ms in enumerate(mu_psi):
        #             arr[i][0:ms.shape[0], 0:ms.shape[1]] = ms
        #         kwargs['mu_psi'] = arr
        #     super().add(**kwargs)

    def heterogen(self, lsv_id):
        """
        Accessor function for PSI file.

        :param lsv_id:
        :return:
        """
        return self._Heterogen(self, lsv_id)
