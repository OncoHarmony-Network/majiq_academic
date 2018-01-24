import numpy

from voila.api.matrix_hdf5 import DeltaPsi, Psi
from voila.vlsv import get_expected_dpsi, VoilaLsv


def unpack_means(value):
    if numpy.size(value, 0) == 1:
        value = numpy.append(value, numpy.array(1 - value[0]))
    return value


def unpack_bins(value):
    if numpy.size(value, 0) == 1:
        value = numpy.append(value, [numpy.flip(value[-1], 0)], axis=0)
    return value


class ViewPsi(Psi):
    class _ViewPsi(Psi._Psi):
        def get(self):
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'coordinates', tuple(self.coordinates)

            for key in self.fields:
                if key == 'bins':
                    yield 'group_bins', dict(self.group_bins)
                elif key == 'means':
                    yield 'group_means_rounded', dict(self.group_means_rounded)
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

            yield from map(int, coords)

        @property
        def means(self):
            return unpack_means(next(super().get('means'))[1])

        @property
        def group_means(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.means

        @property
        def group_means_rounded(self):
            for group_name, means in self.group_means:
                yield group_name, numpy.around(means, decimals=3)

        @property
        def bins(self):
            return unpack_bins(next(super().get('bins'))[1])

        @property
        def group_bins(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.bins

        @property
        def variances(self):
            means = self.means
            for bin in self.bins:
                step_bins = 1.0 / len(bin)
                projection_prod = bin * numpy.arange(step_bins / 2, 1, step_bins) ** 2
                yield numpy.sum(projection_prod) - means[-1] ** 2

        @property
        def lsv_type(self):
            if self._lsv_type is None:
                self._lsv_type = next(super().get('lsv_type'))[1]
            return self._lsv_type

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
        return self._ViewPsi(self, lsv_id)


class ViewDeltaPsi(DeltaPsi):
    class _ViewDeltaPsi(DeltaPsi._DeltaPsi):
        def get(self):
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'coordinates', self.coordinates
            yield 'excl_incl', self.excl_incl

            for key in self.fields:
                if key == 'group_bins':
                    yield key, dict(self.group_bins)
                elif key == 'group_means':
                    yield 'group_means_rounded', dict(self.group_means_rounded)
                elif key == 'bins':
                    yield 'bins', self.bins
                    yield 'means_rounded', self.means_rounded
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
        def bins(self):
            return unpack_bins(next(super().get('bins'))[1])

        @property
        def group_bins(self):
            group_names = self.matrix_hdf5.group_names
            group_bins = next(super().get('group_bins'))[1]
            for group_name, value in zip(group_names, group_bins):
                yield group_name, unpack_bins(value)

        @property
        def means(self):
            for b in self.bins:
                yield get_expected_dpsi(b)

        @property
        def means_rounded(self):
            return numpy.around(tuple(self.means), decimals=3)

        @property
        def group_means(self):
            for value in next(super().get('group_means'))[1]:
                yield unpack_means(value)

        @property
        def group_means_rounded(self):
            group_names = self.matrix_hdf5.group_names
            for group_name, means in zip(group_names, self.group_means):
                yield group_name, numpy.around(means, decimals=3)

        @property
        def excl_incl(self):
            for mean in self.means:
                if mean < 0:
                    yield [-mean, 0]
                else:
                    yield [0, mean]

    def delta_psi(self, lsv_id):
        return self._ViewDeltaPsi(self, lsv_id)


class ViewMatrix(ViewDeltaPsi, ViewPsi):
    def get_gene_ids(self, args=None):
        if args:
            if args.gene_ids:
                yield from args.gene_ids
            elif args.lsv_ids:
                for lsv_id in args.lsv_ids:
                    yield self.lsv_id_to_gene_id(lsv_id)

        yield from self.h['lsvs'].keys()

    def get_lsv_means(self, lsv_id):
        bins = next(self.get(lsv_id, 'bins'))[1]
        for b in bins:
            yield get_expected_dpsi(b)

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
            for lsv_id in self.h['lsvs'][gene_id].keys():
                # Search for LSV
                if not lsv_ids or lsv_id in lsv_ids:
                    if not threshold or VoilaLsv.is_lsv_changing(self.get_lsv_means(lsv_id), threshold):
                        yield lsv_id

    @property
    def metadata(self):
        metadata = super().metadata
        experiment_names = metadata['experiment_names']
        group_names = super().group_names

        metadata['experiment_names'] = [numpy.insert(exps, 0, '{0} Combined'.format(group)) for exps, group in
                                        zip(experiment_names, group_names)]

        return metadata
