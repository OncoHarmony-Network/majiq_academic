import math
import os
from itertools import islice

import h5py
import numpy as np

from voila import constants
from voila.api.splice_graph_type import Gene, Junction, Exon
from voila.api.voila_hdf5 import VoilaHDF5
from voila.utils.exceptions import GeneIdNotFoundInVoilaFile
from voila.utils.voila_log import voila_log


class SpliceGraphHDF5:
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
        except (ValueError, RuntimeError):
            del self.hdf5

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
        pass

    def add_experiment_names(self, experiment_names):
        self.hdf5.attrs['experiment_names'] = np.array(experiment_names, dtype=h5py.special_dtype(vlen=np.unicode))

    def get_experiments(self):
        return self.hdf5.attrs['experiment_names']


class Genes(SpliceGraphHDF5):
    def gene(self, id, **kwargs):
        return Gene(self.hdf5_grp('Gene', id), **kwargs)

    @property
    def genes(self):
        gene_hdf5 = self.hdf5['Gene']
        return (Gene(gene_hdf5[id]) for id in gene_hdf5)

    def combined_genes(self, experiments):
        all_experiments = list(self.get_experiments())
        gene_dict = {}

        for gene in self.genes:

            comb_gene = None

            for exp in experiments:

                exp_idx = all_experiments.index(exp)

                if comb_gene is None:
                    comb_gene = gene.get_experiment(exp_idx)
                else:
                    new_gene = gene.get_experiment(exp_idx)
                    for x, y in zip(comb_gene['junctions'], new_gene['junctions']):
                        x['reads'] += y['reads']

            gene_dict[gene.id] = comb_gene

        return gene_dict

    def get_gene_ids(self, args=None):
        if args and args.gene_ids:
            return args.gene_ids

        if args and hasattr(args, 'lsv_ids') and args.lsv_ids:
            return (lsv_id.split(':')[0] for lsv_id in args.lsv_ids)

        return self.hdf5['Gene'].keys()

    def get_page_count(self, args):
        gene_count = 0
        log = voila_log()

        log.debug('Start page count')

        if hasattr(args, 'voila_file'):
            with VoilaHDF5(args.voila_file, 'r') as v:
                for gene_id in self.get_gene_ids(args):
                    try:
                        if any(v.get_lsvs(args, gene_id)):
                            gene_count += 1
                    except GeneIdNotFoundInVoilaFile:
                        pass

        else:
            log.debug('Gene limit is set to {0}'.format(args.limit))
            for _ in self.get_gene_ids(args):
                gene_count += 1
                if gene_count == args.limit:
                    break

        log.debug('End page count')

        return int(math.ceil(gene_count / float(constants.MAX_GENES)))

    def get_paginated_genes_with_lsvs(self, args):
        log = voila_log()
        log.debug('Getting paginated genes with LSVs')

        gene_list = []
        lsv_dict = {}

        with VoilaHDF5(args.voila_file, 'r') as v:
            for gene_id in self.get_gene_ids(args):
                try:
                    lsvs = tuple(v.get_lsvs(args, gene_id=gene_id))
                except GeneIdNotFoundInVoilaFile:
                    lsvs = None

                if lsvs:
                    lsv_dict[gene_id] = tuple(v.get_voila_lsv(gene_id, lsv_id) for gene_id, lsv_id in lsvs)
                    gene_list.append(gene_id)

                if len(gene_list) == constants.MAX_GENES:
                    yield lsv_dict, gene_list
                    gene_list = []
                    lsv_dict = {}

            if gene_list:
                yield lsv_dict, gene_list

    def get_paginated_genes(self, args):
        log = voila_log()
        log.debug('Getting paginated genes')

        gene_list = []
        gene_count = 0

        for gene_id in self.get_gene_ids(args):
            log.debug('Found {0}'.format(gene_id))
            gene_list.append(gene_id)
            gene_count += 1

            if gene_count == args.limit:
                break

            if len(gene_list) == constants.MAX_GENES:
                yield gene_list
                gene_list = []

        if gene_list:
            yield gene_list


class Junctions(SpliceGraphHDF5):
    def junction(self, id, **kwargs):
        return Junction(self.hdf5_grp('Junctions', id), **kwargs)

    @property
    def junctions(self):
        juncs_hdf5 = self.hdf5['Junctions']
        return (Junction(juncs_hdf5[id]) for id in juncs_hdf5)


class Exons(SpliceGraphHDF5):
    def exon(self, id, **kwargs):
        return Exon(self.hdf5_grp('Exons', id), **kwargs)

    @property
    def exons(self):
        exons_hdf5 = self.hdf5['Exons']
        return (Exon(exons_hdf5[id]) for id in exons_hdf5)
