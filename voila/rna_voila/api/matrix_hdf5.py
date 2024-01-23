import csv
from abc import ABC, abstractmethod
from typing import List

import h5py
import numpy as np

from rna_voila import constants
from rna_voila.api.matrix_utils import generate_excl_incl, generate_means, generate_high_probability_non_changing, \
    generate_variances, generate_standard_deviations
from rna_voila.exceptions import LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile
from rna_voila.vlsv import collapse_matrix, matrix_area


def lsv_id_to_gene_id(lsv_id):
    return ':'.join(lsv_id.split(':')[:-2])

opened_voila_files = {}

def _open_hdf5(filename, mode):
    return h5py.File(filename, mode, libver='latest')

def open_hdf5(filename, mode):
    from rna_voila.config import ViewConfig
    from rna_voila.voila_log import voila_log
    config = ViewConfig()

    filename = str(filename)
    if mode == 'r':
        if config.memory_map_hdf5 or config.preserve_handles_hdf5:
            if filename not in opened_voila_files:
                if config.memory_map_hdf5:
                    voila_log().debug(f'memory mapping {filename}')
                    opened_voila_files[filename] = h5py.File(filename, mode, libver='latest', driver='core', backing_store=False)
                else:
                    voila_log().debug(f'opening handle {filename}')
                    opened_voila_files[filename] = h5py.File(filename, mode, libver='latest', swmr=True)
            return opened_voila_files[filename]

    return _open_hdf5(filename, mode)

class MatrixHdf5:
    LSVS = 'lsvs'

    def __init__(self, filename, mode='r', voila_file=True, voila_tsv=False, pre_config=False):
        """
        Access voila's HDF5 file.

        :param filename: name of voila file
        :param mode: generally r or w
        """
        self.voila_tsv = voila_tsv
        self.voila_file = voila_file
        self.mode = mode
        self.dt = h5py.special_dtype(vlen=str)
        self._group_names = None
        self._tsv_writer = None
        self._tsv_file = None
        self._filename = filename
        self._prior = None

        try:
            from rna_voila.config import ViewConfig
            ViewConfig()
        except KeyError:
            pre_config = True

        self._pre_config = pre_config

        if voila_file:
            self.h = _open_hdf5(filename, mode) if pre_config else open_hdf5(filename, mode)


    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def has_index(self):
        return 'index' in self.h

    def write_index(self, index, hashval):
        if 'index' in self.h:
            del self.h['index']
        if 'input_hash' in self.h:
            del self.h['input_hash']
        self.h.create_dataset('index', index.shape, data=index)
        self.h.create_dataset("input_hash", (1,), dtype="S40", data=(hashval.encode('utf-8'),))

    def get_index(self):
        return self.h['index'][()]

    def close(self):
        from rna_voila.config import ViewConfig
        if self.voila_file:

            if self._pre_config:
                self.h.close()
            elif not ViewConfig().memory_map_hdf5 and not ViewConfig().preserve_handles_hdf5:
                self.h.close()
            elif self.mode != 'r':
                self.h.close()

        if self.voila_tsv:
            self._tsv_file.close()

    def add_dataset(self, *args, **kwargs):
        """
        Add dataset to hdf5 file using h5py api.
        :param args: tree key values
        :param kwargs: kwargs for the h5py create_dataset method.
        :return: None
        """
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
        """
        Add row to tsv file.  This is used in the quantification step of the pipeline.

        :param matrix_type: Matrix class for analysis type.
        :param kwargs: data for voila file
        :return: None
        """
        if self._tsv_writer is None:
            tsv_filename = '.'.join(self._filename.split('.')[:-1]) + '.tsv'
            self._tsv_file = open(tsv_filename, 'w', newline='')
            fieldnames = matrix_type.tsv_fieldnames()
            self._tsv_writer = csv.DictWriter(self._tsv_file, fieldnames=fieldnames, delimiter='\t')
            self._tsv_writer.writeheader()

        junctions = kwargs.get('junctions').tolist()
        lsv_type = kwargs.get('lsv_type', '')
        intron_ret = int(lsv_type[-1] == 'i')

        junc_coords = junctions
        if intron_ret:
            ir_coords = junctions[-1:]
        else:
            ir_coords = []

        row = matrix_type.tsv_row(**kwargs)

        row.update({
            'gene_id': lsv_id_to_gene_id(matrix_type.lsv_id),
            'lsv_id': matrix_type.lsv_id,
            'lsv_type': lsv_type,
            'num_junctions': len(junctions) - intron_ret,
            'num_exons': matrix_type.exon_count,
            'junctions_coords': ';'.join('{0}-{1}'.format(start, end) for start, end in junc_coords),
            'ir_coords': ';'.join('{0}-{1}'.format(start, end) for start, end in ir_coords)
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
            try:
                raise GeneIdNotFoundInVoilaFile(self.h.filename, gene_id)
            except:
                raise GeneIdNotFoundInVoilaFile(None, gene_id)

        try:
            lsv_grp = gene_grp[lsv_id]
        except KeyError:
            try:
                raise LsvIdNotFoundInVoilaFile(self.h.filename, lsv_id)
            except:
                raise LsvIdNotFoundInVoilaFile(None, lsv_id)

        try:
            return lsv_grp[key][()]
        except KeyError:

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
            yield key, lsv_grp[key][()]

    @property
    def prior(self):
        """
        Get Prior from hdf5 file if it hasn't already been cached.
        :return: numpy matrix
        """

        if self._prior is None:
            self._prior = self.h['metadata']['prior'][()]
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
        """
        Gets analysis type from h5py file.
        :return:
        """
        return self.h['metadata']['analysis_type'].asstr()[()]

    @analysis_type.setter
    def analysis_type(self, a):
        """
        Saved analysis type to h5py file.
        :param a: analysis type string
        :return: None
        """
        self.create_dataset('metadata/analysis_type', data=a)

    @property
    def group_names(self):
        """
        If group names haven't already been cached, then cache them before returning list of group names.
        :return: list of strings
        """
        if self._group_names is None:
            self._group_names = self.h['metadata']['group_names'].asstr()[()].tolist()
        return self._group_names

    @group_names.setter
    def group_names(self, n):
        """
        Cache group names and then save them to the h5py file.
        :param n: list of group names
        :return: None
        """
        self._group_names = n
        self.create_dataset('metadata/group_names', data=np.array(n, dtype=self.dt))

    @property
    def splice_graph_experiment_names(self):
        """
        List of experiment names for creating splice graphs.  This will be same list of experiment names as is saved to
        the voila file, but if there are more then one experiment name, there will be a 'Combined' prepended to the list
        of experiment names.
        :return: list of experiment names
        """
        exp_names = []
        for grp, exp in zip(self.group_names, self.experiment_names):
            if len(exp) > 1:
                exp_names.append([grp + ' Combined'] + exp)
            else:
                exp_names.append(exp)
        return exp_names

    @property
    def experiment_names(self):
        """
        Get list of experiment names using h5py api.
        :return:
        """
        return self.h['metadata']['experiment_names'].asstr()[()].tolist()

    @experiment_names.setter
    def experiment_names(self, ns):
        """
        Set list of experiment names to the h5py file. Some care has be taken when the list of experiment names isn't
        even across all columns.  For instance, if one group has more experiments then the rest of the groups.
        :param ns: list of experiment names
        :return: None
        """
        ns = [[e.decode('utf-8') for e in es] for es in ns]
        arr = np.empty((2, max(len(n) for n in ns)), dtype=self.dt)
        arr.fill('')
        for i, n in enumerate(ns):
            arr[i][0:len(n)] = n
        self.create_dataset('metadata/experiment_names', data=arr)

    @property
    def stat_names(self):
        """
        List of stats used in this quantification.
        :return: list of strings
        """
        return self.h['metadata']['stat_names'].asstr()[()]

    @stat_names.setter
    def stat_names(self, s):
        """
        Set list of stat names.
        :param s: list of stat names
        """
        self.h.create_dataset('metadata/stat_names', data=np.array(s, dtype=self.dt))

    @property
    def psi_samples(self):
        """
        Number of psi samples used for p-value quantile computations

        :return: int, number of psisamples used in MAJIQ HET computation
        :raises: KeyError if used before set or reading old VOILA file before
        this field was added
        """
        return self.h['metadata']['psi_samples'][()]

    @psi_samples.setter
    def psi_samples(self, psi_samples):
        """
        Record number of psi samples used for p-value quantile computations

        :param psi_samples: int, parameter used in majiq het for number of psi
        posterior samples used for calculating test statistics
        """
        self.h.create_dataset(
            'metadata/psi_samples', data=np.array(psi_samples, dtype=int)
        )

    @property
    def test_percentile(self):
        """
        Percentile of p-values from psi samples reported for test_quantiles

        :note: it's actually the *quantile*, not the percentile, but is named
        to match nomenclature used in MAJIQ
        :return: float, value used in majiq het for this parameter
        :raises: KeyError if used before set or reading old VOILA file before
        this field was added
        """
        return self.h['metadata']['test_percentile'][()]

    @test_percentile.setter
    def test_percentile(self, test_percentile):
        """
        Record number of psi samples used for p-value quantile computations

        :param psi_samples: float, parameter used in majiq het for test_percentile
        """
        self.h.create_dataset(
            'metadata/test_percentile', data=np.array(test_percentile, dtype=float)
        )

    @property
    def gene_ids(self):
        """
        A generator of all the gene ids in this file.
        :return: generator
        """
        yield from self.h['lsvs']

    @property
    def num_lsv_ids(self):
        return len(self.h['lsvs'])

    def lsv_ids(self, gene_ids=None):
        """
        A generator of lsv ids. If gene_ids is set, then only the lsvs for those gene ids will be provided.
        :param gene_ids: list of gene ids
        :return: generator
        """
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
        """
        File version as set using value from contants.
        :return:
        """
        metadata = self.h['metadata']
        try:
            return metadata['file_version'][()]
        except KeyError:
            return -1

    @file_version.setter
    def file_version(self, version):
        """
        File version set from constants.
        :param version:
        :return:
        """
        self.create_dataset('metadata/file_version', data=version)

    def create_dataset(self, *args, **kwargs):
        """
        Wrapper function for creating datasets with h5py api.
        :param args: arguments
        :param kwargs: keyword arguments
        :return: None
        """
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

    def get(self, key):
        """
        Retrieve value from Matrix file.

        :param key: Key for value
        :return: Value
        """
        return self.matrix_hdf5.get(self.lsv_id, key)

    def get_many(self, keys):
        """
        Retrieve many values using list of keys

        :param keys: list of keys for values
        :return: Generator of key, values
        """
        return self.matrix_hdf5.get_many(self.lsv_id, keys)

    @property
    def exists(self):
        """
        Does lsv id exist in h5py file.
        :return: boolean
        """
        gene_id = self.gene_id
        lsv_id = self.lsv_id
        h = self.matrix_hdf5.h
        return gene_id in h['lsvs'] and lsv_id in h[gene_id]

    @property
    def lsv_type(self):
        """
        Get lsv type from h5py file.
        :return: string
        """
        if self._lsv_type is None:
            self._lsv_type = self.get('lsv_type').decode()
        return self._lsv_type

    @property
    def intron_retention(self):
        """
        Using the lsv type, does this lsv have intron retention.
        :return: boolean
        """
        return 'i' == self.lsv_type[-1]

    @property
    def reference_exon(self):
        """
        Get coordinates for the lsv reference exon.
        :return:
        """
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
        return tuple(ref_exon)

    @property
    def target(self):
        """
        Use lsv_id to check if lsv is target.
        :return: boolean
        """
        return self.lsv_id.split(':')[-2] == 't'

    @property
    def source(self):
        """
        Inverse of target.
        :return: boolean
        """
        return not self.target

    @property
    def gene_id(self):
        """
        Convert lsv_id to gene id.
        :return: string
        """
        return lsv_id_to_gene_id(self.lsv_id)

    @property
    def a5ss(self):
        """
        Using lsv type, does this lsv have 5 prime splice sites.
        :return: boolean
        """
        return 'A5SS' in [self.reference_exon_ss(), self.other_exons_ss()]

    @property
    def a3ss(self):
        """
        Using lsv type, does this lsv have 3 prime splice sites.
        :return: boolean
        """
        return 'A3SS' in [self.reference_exon_ss(), self.other_exons_ss()]

    def reference_exon_ss(self):
        """
        Check for 3 prime or 5 prime splice sites in reference exon.
        :return: list of strings
        """
        try:

            ss = filter(lambda x: x != 'i', self.lsv_type.split('|')[1:])
            ss = map(lambda x: x.split('.')[0].split('e')[0], ss)

            if len(set(ss)) > 1:
                if self.lsv_type[0] == 's':
                    return 'A5SS'
                else:
                    return 'A3SS'

        except IndexError:

            if constants.NA_LSV in self.lsv_type:
                return constants.NA_LSV
            raise

    def other_exons_ss(self):
        """
        Find 3 prime or 5 prime splice sites in exons that aren't the reference exon.
        :return: List of strings
        """
        try:

            ss = filter(lambda lt: lt != 'i', self.lsv_type.split('|')[1:])
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
            if constants.NA_LSV in self.lsv_type:
                return constants.NA_LSV
            raise

    @property
    def exon_skipping(self):
        """
        Using lsv type, does this lsv have exon skipping.
        :return: boolean
        """
        try:
            return self.exon_count > 2
        except TypeError:
            if self.exon_count == constants.NA_LSV:
                return constants.NA_LSV
            raise

    @property
    def exon_count(self):
        """
        Using lsv type, how many exons are in this lsv.
        :return: integer
        """
        try:
            exons = filter(lambda x: x != 'i', self.lsv_type.split('|')[1:])
            exons = map(lambda x: x.split('.')[0].split('e')[1], exons)
            return len(set(exons)) + 1
        except IndexError:
            if constants.NA_LSV in self.lsv_type:
                return constants.NA_LSV
            raise

    @property
    def binary(self):
        """
        Using lsv type, is this lsv binary.
        :return: boolean
        """
        return len(self.lsv_type.split('|')[1:]) == 2

    @property
    def complex(self):
        """
        Inverse of binary.
        :return: boolean
        """
        return not self.binary

    @property
    def junctions(self):
        """
        Including intron retention, a 2d list of junction coordinates. Intron retention is the last set of coordinates.
        :return: numpy matrix
        """
        return self.get('junctions')


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
            """
            Delta PSI TSV fieldnames.
            :return: list of fieldnames
            """
            group_names = self.matrix_hdf5.group_names
            return [
                'gene_id',
                'lsv_id',
                'lsv_type',
                'mean_dpsi_per_lsv_junction',
                'probability_changing',
                'probability_non_changing',
                '{}_mean_psi'.format(group_names[0]),
                '{}_mean_psi'.format(group_names[1]),
                'num_junctions',
                'num_exons',
                'junctions_coords',
                'ir_coords'
            ]

        def tsv_row(self, **kwargs):
            """
            Add columns that are not general to all analysis types.
            :param kwargs: voila file data
            :return: dictionary
            """
            bins = kwargs['bins']
            means = generate_means(bins)
            excl_incl = generate_excl_incl(means)
            threshold = 0.2
            non_changing_threshold = 0.05

            row = {
                'mean_dpsi_per_lsv_junction': ';'.join(
                    f'{excl_incl[i][1] - excl_incl[i][0]:0.4f}' for i in range(np.size(bins, 0))),
                'probability_changing': ';'.join(f'{matrix_area(b, threshold):0.3e}' for b in bins),
                'probability_non_changing': ';'.join(
                    map(
                        (lambda x: f'{x:0.3e}'),
                        generate_high_probability_non_changing(
                            self.intron_retention,
                            self.matrix_hdf5.prior,
                            non_changing_threshold, bins
                        )
                    )
                ),
            }

            for group_name, means in zip(self.matrix_hdf5.group_names, kwargs['group_means']):
                row[f'{group_name}_mean_psi'] = ';'.join('%0.4f' % i for i in means)

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

        @staticmethod
        def tsv_fieldnames():
            """
            PSI fieldnames.
            :return: list of psi fieldnames
            """
            return [
                'gene_id',
                'lsv_id',
                'lsv_type',
                'mean_psi_per_lsv_junction',
                'stdev_psi_per_lsv_junction',
                'num_junctions',
                'num_exons',
                'junctions_coords',
                'ir_coords',
            ]

        @staticmethod
        def tsv_row(**kwargs):
            """
            Add columns that are not general to all analysis types.
            :param kwargs: voila file data
            :return:  dictionary
            """
            bins = kwargs['bins']
            return {
                'mean_psi_per_lsv_junction': ';'.join(map((lambda x: f'{x:0.4f}'), kwargs['means'])),
                'stdev_psi_per_lsv_junction': ';'.join(map((lambda x: f'{x:0.4f}'), generate_standard_deviations(bins)))
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
