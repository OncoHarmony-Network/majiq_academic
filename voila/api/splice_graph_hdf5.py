import os
from itertools import islice

import h5py

from voila.api.splice_graph_type import Gene, Junction, Exon


class SpliceGraphHDF5():
    def __init__(self, filename, mode='r'):
        self.hdf5 = None

        try:
            self.hdf5 = h5py.File(filename, mode=mode, swmr=True)
        except ValueError as ve:
            if str(ve) == 'The SWMR feature is not available in this version of the HDF5 library':
                self.hdf5 = h5py.File(filename, mode=mode)
            else:
                raise

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        pass

    def close(self):
        try:
            self.hdf5.flush()
            self.hdf5.close()
        except ValueError:
            pass

    def hdf5_grp(self, t, id):
        id = str(id)
        try:
            return self.hdf5[t][id]
        except KeyError:
            return self.hdf5.create_group(os.path.join(t, id))

    @staticmethod
    def get(itr, n):
        try:
            return next(islice(itr, n, n + 1))
        except StopIteration:
            return None

    def check_version(self):
        # todo: check file version
        return ''

    def add_experiment_names(self, experiment_names):
        # todo: add experiment names to metainfo
        pass


class Genes(SpliceGraphHDF5):
    def gene(self, id, **kwargs):
        return Gene(self.hdf5_grp('Gene', id), **kwargs)

    def genes(self):
        gene_hdf5 = self.hdf5['Gene']
        return (Gene(gene_hdf5[id]) for id in gene_hdf5)


class Junctions(SpliceGraphHDF5):
    def junction(self, id, **kwargs):
        return Junction(self.hdf5_grp('Junctions', id), **kwargs)

    def junctions(self):
        juncs_hdf5 = self.hdf5['Junctions']
        return (Junction(juncs_hdf5[id]) for id in juncs_hdf5)


class Exons(SpliceGraphHDF5):
    def exon(self, id, **kwargs):
        return Exon(self.hdf5_grp('Exons', id), **kwargs)

    def exons(self):
        exons_hdf5 = self.hdf5['Exons']
        return (Exon(exons_hdf5[id]) for id in exons_hdf5)
