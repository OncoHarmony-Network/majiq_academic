import os
from abc import ABC, abstractmethod

import h5py
import numpy

from voila.vlsv import collapse_matrix


class MatrixHdf5:
    def __init__(self, filename, mode='r'):
        self.h = h5py.File(filename, mode, libver='latest')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.h.close()

    def add(self, lsv_id, key, value):
        gene_id = self.lsv_id_to_gene_id(lsv_id)
        self.h.create_dataset(os.path.join('lsvs', gene_id, lsv_id, key), data=value)

    def add_multiple(self, lsv_id, **kwargs):
        for key, value in kwargs.items():
            self.add(lsv_id, key, value)

    def get(self, lsv_id, *args):
        gene_id = self.lsv_id_to_gene_id(lsv_id)
        lsv_grp = self.h['lsvs'][gene_id][lsv_id]

        if args:
            keys = args
        else:
            keys = lsv_grp.keys()

        for key in keys:
            yield key, lsv_grp[key].value

    @staticmethod
    def lsv_id_to_gene_id(lsv_id):
        return ':'.join(lsv_id.split(':')[:-2])

    @property
    def analysis_type(self):
        return self.h['metadata']['analysis_type'].value

    @analysis_type.setter
    def analysis_type(self, analysis_type):
        self.h.create_dataset('metadata/analysis_type', data=analysis_type)

    @property
    def group_names(self):
        return self.h['metadata']['group_names'].value

    @group_names.setter
    def group_names(self, names):
        dt = h5py.special_dtype(vlen=numpy.unicode)
        self.h.create_dataset(os.path.join('metadata', 'group_names'),
                              data=numpy.array(names, dtype=dt))

    @property
    def experiment_names(self):
        return self.h['metadata']['experiment_names'].value

    @experiment_names.setter
    def experiment_names(self, names):
        dt = h5py.special_dtype(vlen=numpy.unicode)
        self.h.create_dataset(os.path.join('metadata', 'experiment_names'),
                              data=numpy.array(names, dtype=dt))

    @property
    def metadata(self):
        return {
            'analysis_type': self.analysis_type,
            'group_names': self.group_names,
            'experiment_names': self.experiment_names
        }


class MatrixType(ABC):
    @abstractmethod
    def __init__(self, matrix_hdf5, lsv_id, fields):
        self.matrix_hdf5 = matrix_hdf5
        self.lsv_id = lsv_id
        self.fields = fields
        self._lsv_type = None

    def add(self, **kwargs):
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

    def get(self, *args):
        return self.matrix_hdf5.get(self.lsv_id, *args)


class DeltaPsi(MatrixHdf5):
    class _DeltaPsi(MatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            fields = ('bins', 'group_bins', 'group_means', 'lsv_type')
            super().__init__(matrix_hdf5, lsv_id, fields)

        def add(self, **kwargs):
            bins_list = kwargs.get('bins', [])
            bins = [collapse_matrix(bins) for bins in bins_list]
            kwargs['bins'] = bins
            super().add(**kwargs)

    def delta_psi(self, lsv_id):
        return self._DeltaPsi(self, lsv_id)


class Psi(MatrixHdf5):
    class _Psi(MatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            fields = ('bins', 'means', 'lsv_type')
            super().__init__(matrix_hdf5, lsv_id, fields)

    def psi(self, lsv_id):
        return self._Psi(self, lsv_id)
