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
        def get(self, *args):
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'coordinates', tuple(self.coordinates)
            yield 'gene_id', self.gene_id

            if not args:
                args = self.fields

            for key in args:
                if key == 'bins':
                    yield 'group_bins', dict(self.group_bins)
                elif key == 'means':
                    yield 'group_means_rounded', dict(self.group_means_rounded)
                else:
                    yield from super().get(key)

        @property
        def means(self):
            yield from unpack_means(next(super().get('means'))[1])

        @property
        def group_means(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.means

        @property
        def group_means_rounded(self):
            for group_name, means in self.group_means:
                yield group_name, numpy.around(tuple(means), decimals=3)

        @property
        def bins(self):
            return unpack_bins(next(super().get('bins'))[1])

        @property
        def group_bins(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.bins

        @property
        def variances(self):
            means = tuple(self.means)
            for idx, b in enumerate(self.bins):
                step_bins = 1.0 / b.size
                projection_prod = b * numpy.arange(step_bins / 2, 1, step_bins) ** 2
                yield numpy.sum(projection_prod) - means[idx] ** 2
            # variances = []
            # means = tuple(self.means)
            # for idx, bin in enumerate(self.bins):
            #     step_bins = 1.0 / bin.size
            #     projection_prod = bin * numpy.arange(step_bins / 2, 1, step_bins) ** 2
            #     variances.append(numpy.sum(projection_prod) - (means[idx] ** 2))
            # return variances

        @property
        def junction_count(self):
            return len(tuple(self.means))

    def psi(self, lsv_id):
        return self._ViewPsi(self, lsv_id)


class ViewDeltaPsi(DeltaPsi):
    class _ViewDeltaPsi(DeltaPsi._DeltaPsi):
        def get(self, *args):
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'coordinates', self.coordinates
            yield 'excl_incl', self.excl_incl

            if not args:
                args = self.fields

            for key in args:
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

        @property
        def junction_count(self):
            return len(tuple(self.means))

    def delta_psi(self, lsv_id):
        return self._ViewDeltaPsi(self, lsv_id)


class ViewMatrix(ViewDeltaPsi, ViewPsi):

    def get_lsv_count(self, args):
        return sum(1 for _ in self.get_lsvs(args))

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
            for lsv_id in self.lsv_ids(gene_id):
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
