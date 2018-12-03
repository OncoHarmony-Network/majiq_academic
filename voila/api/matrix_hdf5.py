import csv
from abc import ABC, abstractmethod
from typing import List

import h5py
import numpy as np

from voila import constants
from voila.api.matrix_utils import generate_excl_incl, generate_means, generate_high_probability_non_changing, \
    generate_variances
from voila.exceptions import LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile
from voila.vlsv import collapse_matrix, matrix_area


def lsv_id_to_gene_id(lsv_id):
    return ':'.join(lsv_id.split(':')[:-2])


class MatrixHdf5:
    LSVS = 'lsvs'

    def __init__(self, filename, mode='r', voila_file=True, voila_tsv=False):
        """
        Access voila's HDF5 file.

        :param filename: name of voila file
        :param mode: generally r or w
        """
        self.voila_tsv = voila_tsv
        self.voila_file = voila_file
        self.dt = h5py.special_dtype(vlen=np.unicode)
        self._group_names = None
        self._tsv_writer = None
        self._tsv_file = None
        self._filename = filename
        self._prior = None

        if voila_file:
            self.h = h5py.File(filename, mode, libver='latest')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.voila_file:
            self.h.close()

        if self.voila_tsv:
            self._tsv_file.close()

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

    def add_tsv_row(self, matrix_type, **kwargs):

        if self._tsv_writer is None:
            tsv_filename = '.'.join(self._filename.split('.')[:-1]) + '.tsv'
            self._tsv_file = open(tsv_filename, 'w', newline='')
            fieldnames = matrix_type.tsv_fieldnames()
            self._tsv_writer = csv.DictWriter(self._tsv_file, fieldnames=fieldnames, delimiter='\t')
            self._tsv_writer.writeheader()

        junctions = kwargs.get('junctions').tolist()
        lsv_type = kwargs.get('lsv_type', '')
        intron_ret = int(lsv_type[-1] == 'i')

        if intron_ret:
            junc_coords = junctions[:-1]
            ir_coords = junctions[-2:-1]
        else:
            junc_coords = junctions
            ir_coords = []

        row = matrix_type.tsv_row(**kwargs)

        row.update({
            'Gene ID': lsv_id_to_gene_id(matrix_type.lsv_id),
            'LSV ID': matrix_type.lsv_id,
            'LSV Type': lsv_type,
            'A5SS': matrix_type.a5ss,
            'A3SS': matrix_type.a3ss,
            'ES': matrix_type.exon_skipping,
            'Num. Junctions': len(junctions) - intron_ret,
            'Num. Exons': matrix_type.exon_count,
            'Junctions coords': ';'.join('{0}-{1}'.format(start, end) for start, end in junc_coords),
            'IR coords': ';'.join('{0}-{1}'.format(start, end) for start, end in ir_coords)
        })

        self._tsv_writer.writerow(row)

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
        except (KeyError, ValueError):
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
        if self._prior is None:
            self._prior = self.h['metadata']['prior'].value
        return self._prior

    @prior.setter
    def prior(self, ps):
        """
        Assumed that this is received in log space.
        :param ps:
        :return:
        """
        self._prior = [collapse_matrix(p) for p in ps]
        self.create_dataset('metadata/prior', data=self._prior)

    @property
    def analysis_type(self):
        return self.h['metadata']['analysis_type'].value

    @analysis_type.setter
    def analysis_type(self, a):
        self.create_dataset('metadata/analysis_type', data=a)

    @property
    def group_names(self):
        if self._group_names is None:
            self._group_names = self.h['metadata']['group_names'].value.tolist()
        return self._group_names

    @group_names.setter
    def group_names(self, n):
        self._group_names = n
        self.create_dataset('metadata/group_names', data=np.array(n, dtype=self.dt))

    @property
    def splice_graph_experiment_names(self):
        exp_names = []
        for grp, exp in zip(self.group_names, self.experiment_names):
            if len(exp) > 1:
                exp_names.append([grp + ' Combined'] + exp)
            else:
                exp_names.append(exp)
        return exp_names

    @property
    def experiment_names(self):
        return self.h['metadata']['experiment_names'].value.tolist()

    @experiment_names.setter
    def experiment_names(self, ns):
        ns = [[e.decode('utf-8') for e in es] for es in ns]
        arr = np.empty((2, max(len(n) for n in ns)), dtype=self.dt)
        arr.fill('')
        for i, n in enumerate(ns):
            arr[i][0:len(n)] = n
        self.create_dataset('metadata/experiment_names', data=arr)

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
            gene_ids = self.gene_ids

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
        self.create_dataset('metadata/file_version', data=version)

    def create_dataset(self, *args, **kwargs):
        if self.voila_file:
            self.h.create_dataset(*args, **kwargs)


class MatrixType(ABC):
    @abstractmethod
    def __init__(self, matrix_hdf5, lsv_id, fields):
        """
        The analysis type of data found in matrix file.

        :param matrix_hdf5: matrix HDF5 object
        :param lsv_id: unique LSV identifier
        :param fields: fields allowed to be stored in this file
        """
        self.voila_tsv = matrix_hdf5.voila_tsv
        self.voila_file = matrix_hdf5.voila_file
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
        if self.voila_file:
            self.matrix_hdf5.add_multiple(self.lsv_id, **kwargs)

        if self.voila_tsv:
            if self._lsv_type is None:
                self._lsv_type = kwargs.get('lsv_type', None)
            self.matrix_hdf5.add_tsv_row(self, **kwargs)

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
    def exists(self):
        gene_id = self.gene_id
        lsv_id = self.lsv_id
        h = self.matrix_hdf5.h
        return gene_id in h['lsvs'] and lsv_id in h[gene_id]

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
        ref_exon = []
        if len(coords) == 2:
            for coord in coords:
                if coord == 'na':
                    ref_exon.append(-1)
                else:
                    ref_exon.append(int(coord))
        else:
            for coord in coords:
                if coord:
                    if coord == '1':
                        ref_exon.append(-1)
                    else:
                        ref_exon.append(int(coord))
        return ref_exon

    @property
    def target(self):
        return self.lsv_id.split(':')[-2] == 't'

    @property
    def source(self):
        return not self.target

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

    @property
    def complex(self):
        return not self.binary


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

        def tsv_fieldnames(self):
            group_names = self.matrix_hdf5.group_names
            return ['Gene ID', 'LSV ID', 'LSV Type', 'E(dPSI) per LSV junction', 'P(|dPSI|>=0.20) per LSV junction',
                    'P(|dPSI|<=0.05) per LSV junction', '{} E(PSI)'.format(group_names[0]),
                    '{} E(PSI)'.format(group_names[1]), 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons',
                    'Junctions coords', 'IR coords']

        def tsv_row(self, **kwargs):
            bins = kwargs['bins']
            means = generate_means(bins)
            excl_incl = generate_excl_incl(means)
            threshold = 0.2
            non_changing_threshold = 0.05

            row = {
                'E(dPSI) per LSV junction': ';'.join(
                    str(excl_incl[i][1] - excl_incl[i][0]) for i in range(np.size(bins, 0))),
                'P(|dPSI|>=%.2f) per LSV junction' % threshold: ';'.join(str(matrix_area(b, threshold)) for b in bins),
                'P(|dPSI|<=%.2f) per LSV junction' % non_changing_threshold: ';'.join(
                    map(str, generate_high_probability_non_changing(self.intron_retention, self.matrix_hdf5.prior,
                                                                    non_changing_threshold, bins))),
            }

            for group_name, means in zip(self.matrix_hdf5.group_names, kwargs['group_means']):
                row[group_name + ' E(PSI)'] = ';'.join('%.3f' % i for i in means)

            return row

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

        def tsv_fieldnames(self):
            return ['Gene ID', 'LSV ID', 'LSV Type', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction', 'A5SS',
                    'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'Junctions coords', 'IR coords']

        def tsv_row(self, **kwargs):
            bins = kwargs['bins']
            return {
                'E(PSI) per LSV junction': ';'.join(map(str, kwargs['means'])),
                'Var(E(PSI)) per LSV junction': ';'.join(map(str, generate_variances(bins)))
            }

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

    def heterogen(self, lsv_id):
        """
        Accessor function for PSI file.

        :param lsv_id:
        :return:
        """
        return self._Heterogen(self, lsv_id)
