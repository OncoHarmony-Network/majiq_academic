import os
from abc import ABC, abstractmethod

import h5py
import numpy

from voila.vlsv import collapse_matrix


def lsv_id_to_gene_id(lsv_id):
    return ':'.join(lsv_id.split(':')[:-2])


class MatrixHdf5:
    def __init__(self, filename, mode='r'):
        filename = os.path.expanduser(filename)
        self.h = h5py.File(filename, mode, libver='latest')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.h.close()

    def add(self, lsv_id, key, value):
        gene_id = lsv_id_to_gene_id(lsv_id)
        self.h.create_dataset(os.path.join('lsvs', gene_id, lsv_id, key), data=value)

    def add_multiple(self, lsv_id, **kwargs):
        for key, value in kwargs.items():
            self.add(lsv_id, key, value)

    def get(self, lsv_id, key):
        gene_id = lsv_id_to_gene_id(lsv_id)
        lsv_grp = self.h['lsvs'][gene_id][lsv_id]
        return lsv_grp[key].value

    def get_many(self, lsv_id, keys):
        gene_id = lsv_id_to_gene_id(lsv_id)
        lsv_grp = self.h['lsvs'][gene_id][lsv_id]
        for key in keys:
            yield key, lsv_grp[key].value

    @property
    def prior(self):
        return self.h['prior'].value

    @prior.setter
    def prior(self, p):
        self.h.create_dataset('prior', data=collapse_matrix(p))

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
        dt = h5py.special_dtype(vlen=numpy.unicode)
        self.h.create_dataset('metadata/group_names', data=numpy.array(n, dtype=dt))

    @property
    def experiment_names(self):
        return self.h['metadata']['experiment_names'].value

    @experiment_names.setter
    def experiment_names(self, n):
        dt = h5py.special_dtype(vlen=numpy.unicode)
        self.h.create_dataset('metadata/experiment_names', data=numpy.array(n, dtype=dt))

    @property
    def metadata(self):
        return {
            'analysis_type': self.analysis_type,
            'group_names': self.group_names,
            'experiment_names': self.experiment_names
        }

    @property
    def gene_ids(self):
        yield from self.h['lsvs'].keys()

    def lsv_ids(self, gene_id):
        lsvs = self.h['lsvs']
        try:
            yield from lsvs[gene_id].keys()
        except (KeyError, ValueError):
            return ()


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

    def get(self, key):
        return self.matrix_hdf5.get(self.lsv_id, key)

    def get_many(self, keys):
        return self.matrix_hdf5.get_many(self.lsv_id, keys)

    @property
    def lsv_type(self):
        if self._lsv_type is None:
            self._lsv_type = self.get('lsv_type')
        return self._lsv_type

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
        lsv_type = self.lsv_type[2:]
        if lsv_type[-1] == 'i':
            lsv_type = lsv_type[:-2]

        splice_sites = set(j[0] for j in lsv_type.split('|'))
        return len(splice_sites) > 1

    @property
    def _prime3(self):
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
        return self.exon_count > 2

    @property
    def exon_count(self):
        lsv_type = self.lsv_type[2:]
        if lsv_type[-1] == 'i':
            lsv_type = lsv_type[:-2]

        return len(set(x[1:3] for x in lsv_type.split('|'))) + 1


class DeltaPsi(MatrixHdf5):
    class _DeltaPsi(MatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            fields = ('bins', 'group_bins', 'group_means', 'lsv_type', 'junctions')
            super().__init__(matrix_hdf5, lsv_id, fields)

        def add(self, **kwargs):
            bins_list = kwargs.get('bins')
            bins = [collapse_matrix(bins) for bins in bins_list]
            kwargs['bins'] = bins
            junctions = kwargs.get('junctions')
            kwargs['junctions'] = junctions.astype(int)
            super().add(**kwargs)

    def delta_psi(self, lsv_id):
        return self._DeltaPsi(self, lsv_id)


class Psi(MatrixHdf5):
    class _Psi(MatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            fields = ('bins', 'means', 'lsv_type', 'junctions')
            super().__init__(matrix_hdf5, lsv_id, fields)

    def psi(self, lsv_id):
        return self._Psi(self, lsv_id)
