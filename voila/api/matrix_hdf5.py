import os
from abc import ABC, abstractmethod

import h5py
import numpy

from voila.vlsv import VoilaLsv, get_expected_dpsi


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
    def genome(self):
        return self.h['metadata']['genome'].value

    @genome.setter
    def genome(self, genome):
        self.h.create_dataset('metadata/genome', data=genome)

    @property
    def analysis_type(self):
        return self.h['metadata']['analysis_type'].value

    @analysis_type.setter
    def analysis_type(self, analysis_type):
        self.h.create_dataset('metadata/analysis_type', data=analysis_type)

    def add_experiments(self, group_name, experiment_names):
        dt = h5py.special_dtype(vlen=numpy.unicode)
        self.h.create_dataset(os.path.join('metadata', 'experiments', group_name),
                              data=numpy.array([numpy.unicode(e) for e in experiment_names], dtype=dt))

    @property
    def experiments(self):
        experiments = self.h['metadata']['experiments']
        for key in experiments:
            yield key, experiments[key].value

    @property
    def metadata(self):
        return {'genome': self.genome, 'analysis_type': self.analysis_type, 'experiments': dict(self.experiments)}

    def get_gene_ids(self, args=None):
        if args:
            if args.gene_ids:
                return args.gene_ids
            elif args.lsv_ids:
                return (self.lsv_id_to_gene_id(lsv_id) for lsv_id in args.lsv_ids)

        return self.h['lsvs'].keys()

    def get_lsv_ids(self, gene_id):
        yield from self.h['lsvs'][gene_id].keys()

    def get_lsv_means(self, lsv_id):
        bins = dict(self.get(lsv_id, 'bins_list'))['bins_list']
        return tuple(get_expected_dpsi(b) for b in bins)

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
                    try:
                        if not threshold or VoilaLsv.is_lsv_changing(self.get_lsv_means(lsv_id), threshold):
                            yield lsv_id
                    except Exception:
                        print(lsv_id)
                        print(gene_id)
                        raise


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
        def __init__(self, matrix_hdf5, lsv_id, fields):
            super().__init__(matrix_hdf5, lsv_id, fields)

    def delta_psi(self, lsv_id):
        return self._DeltaPsi(self, lsv_id, ('bins', 'psi1', 'psi2', 'means_psi1', 'means_psi2'))


class Psi(MatrixHdf5):
    class _Psi(MatrixType):
        def __init__(self, matrix_hdf5, lsv_id, fields):
            super().__init__(matrix_hdf5, lsv_id, fields)

        def get(self, *args):
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'coordinates', self.coordinates

            if not args:
                args = self.fields

            for key in args:
                if key == 'bins':
                    yield key, self.bins
                elif key == 'means':
                    means = self.means
                    yield key, means
                    yield 'means_rounded', numpy.around(means, decimals=3)
                else:
                    yield from super().get(key)

        @property
        def coordinates(self):
            lsv_coords = self.lsv_id.split(':')[-1]
            if lsv_coords.startswith('-1'):
                coords = (-1, lsv_coords.split('-')[-1])
            elif lsv_coords.endswith('-1'):
                coords = (lsv_coords.split('-')[0], -1)
            else:
                coords = lsv_coords.split('-')

            return list(map(int, coords))

        @property
        def gene_id(self):
            return ':'.join(self.lsv_id.split(':')[:-2])

        @property
        def lsv_type(self):
            if self._lsv_type is None:
                self._lsv_type = next(super().get('lsv_type'))[1]
            return self._lsv_type

        @property
        def means(self):
            value = next(super().get('means'))[1]
            if numpy.size(value, 0) == 1:
                value = numpy.append(value, numpy.array(1 - value[0]))
            return value

        @property
        def bins(self):
            value = next(super().get('bins'))[1]
            if numpy.size(value, 0) == 1:
                value = numpy.append(value, [numpy.flip(value[-1], 0)], axis=0)
            return value

        @property
        def variances(self):
            variances = []
            bins = self.bins
            means = self.means

            if bins is not None and bins.size > 0:
                for bin in bins:
                    step_bins = 1.0 / len(bin)
                    projection_prod = bin * numpy.arange(step_bins / 2, 1, step_bins) ** 2
                    variances.append(numpy.sum(projection_prod) - means[-1] ** 2)
            return variances

        @property
        def target(self):
            return self.lsv_type[0] == 't'

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
        def exon_skipping(self):
            return self.exon_count > 2

        @property
        def exon_count(self):
            lsv_type = self.lsv_type[2:]
            if lsv_type[-1] == 'i':
                lsv_type = lsv_type[:-2]

            return len(set(x[1:3] for x in lsv_type.split('|'))) + 1

        @property
        def junction_count(self):
            return numpy.size(self.means, 0)

    def psi(self, lsv_id):
        return self._Psi(self, lsv_id, ('bins', 'means', 'lsv_type'))

# class MatrixWriterAsync:
#     def __init__(self, filename, mode):
#         self.filename = filename
#         self.mode = mode
#         self.q = Queue()
#         self.e = Event()
#         self.p = Process(target=self._worker, args=(self.q, self.e,))
#         self.p.start()
#
#     def __enter__(self):
#         return self
#
#     def __exit__(self, exc_type, exc_val, exc_tb):
#         self.close()
#
#     def _worker(self, queue, event):
#         with Matrix(self.filename, self.mode) as m:
#             while not (event.is_set() and queue.empty()):
#                 try:
#                     key, value = queue.get_nowait()
#                     m.add(key, value)
#                 except Empty:
#                     pass
#
#     def close(self):
#         self.q.close()
#         self.e.set()
#         self.p.join()
#
#     def add(self, key, value):
#         args = (key, value)
#         try:
#             self.q.put_nowait(args)
#         except Full:
#             self.q.put(args)
