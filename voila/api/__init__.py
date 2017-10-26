from voila.api.splice_graph_sql import Exons, Junctions, Genes
from voila.api.voila_hdf5 import VoilaHDF5


# class OldSpliceGraphs(object):
#     GENES = '/genes'
#     ROOT = '/'
#     VERSION = '/splice_graph_file_version'
#
#     def __init__(self, splice_graph_file_name, mode):
#         """
#         Class for creating and accessing the splice graph file.
#         :param splice_graph_file_name: path to splice graph file
#         :param mode: mode to pass to hdf5
#         """
#         super(OldSpliceGraphs, self).__init__()
#         self.file_name = splice_graph_file_name
#         self.mode = mode
#         self.hdf5 = None
#         self.limit = None
#         self.gene_ids = None
#         self.file_version = None
#
#     def __enter__(self):
#         """
#         Open hdf5 in with block.
#         :return: self
#         """
#         self.hdf5 = h5py.File(self.file_name, self.mode)
#
#         if self.VERSION not in self.hdf5:
#             if self.mode == constants.FILE_MODE.write:
#                 self.hdf5[self.VERSION] = constants.SPLICE_GRAPH_FILE_VERSION
#
#         try:
#             self.file_version = self.hdf5[self.VERSION].value
#         except KeyError:
#             pass
#
#         return self
#
#     def __exit__(self, type, value, traceback):
#         """
#         Close when with block exits.
#         :param type: unused
#         :param value: unused
#         :param traceback: unused
#         :return: None
#         """
#         self.close()
#
#     def close(self):
#         """
#         Close hdf5 file.
#         :return: None
#         """
#         try:
#             self.hdf5.close()
#         except Exception:
#             pass
#
#     def erase_splice_graph_file(self):
#         """
#         Remove splice graph file and reopen it.
#         :return:
#         """
#         os.remove(self.file_name)
#         self.__enter__()
#
#     def add_gene(self, gene):
#         """
#         Add gene object to splice graph file.
#         :param gene: GeneGraphic object
#         :return: None
#         """
#         gene.to_hdf5(self.hdf5)
#
#     def get_page_count(self, args):
#         gene_count = 0
#         log = voila_log()
#
#         log.debug('Start page count')
#
#         if hasattr(args, 'voila_file'):
#             with Voila(args.voila_file, 'r') as v:
#                 for gene_id in self.get_gene_ids(args):
#                     try:
#                         if any(v.get_lsvs(args, gene_id)):
#                             gene_count += 1
#                     except GeneIdNotFoundInVoilaFile:
#                         pass
#
#         else:
#             log.debug('Gene limit is set to {0}'.format(args.limit))
#             for _ in self.get_gene_ids(args):
#                 gene_count += 1
#                 if gene_count == args.limit:
#                     break
#
#         log.debug('End page count')
#
#         return int(math.ceil(gene_count / float(constants.MAX_GENES)))
#
#     def get_gene_name(self, gene_id):
#         return self.hdf5[self.GENES][gene_id].attrs['name']
#
#     def get_gene_ids(self, args=None):
#         if args and args.gene_ids:
#             return args.gene_ids
#
#         if args and hasattr(args, 'lsv_ids') and args.lsv_ids:
#             return (lsv_id.split(':')[0] for lsv_id in args.lsv_ids)
#
#         return self.hdf5[self.GENES].keys()
#
#     def get_genes(self):
#         for gene_id in self.get_gene_ids():
#             yield self.get_gene(gene_id)
#
#     def get_paginated_genes(self, args):
#         log = voila_log()
#         log.debug('Getting paginated genes')
#
#         gene_list = []
#         gene_count = 0
#
#         for gene_id in self.get_gene_ids(args):
#             log.debug('Found {0}'.format(gene_id))
#             gene_list.append(self.get_gene(gene_id))
#             gene_count += 1
#
#             if gene_count == args.limit:
#                 break
#
#             if len(gene_list) == constants.MAX_GENES:
#                 yield gene_list
#                 gene_list = []
#
#         if gene_list:
#             yield gene_list
#
#     def get_paginated_genes_with_lsvs(self, args):
#         log = voila_log()
#         log.debug('Getting paginated genes with LSVs')
#
#         gene_list = []
#         lsv_dict = {}
#
#         with Voila(args.voila_file, 'r') as v:
#             for gene_id in self.get_gene_ids(args):
#                 try:
#                     lsvs = tuple(v.get_lsvs(args, gene_id=gene_id))
#                 except GeneIdNotFoundInVoilaFile:
#                     lsvs = None
#
#                 if lsvs:
#                     gene = self.get_gene(gene_id)
#                     lsv_dict[gene_id] = tuple(v.get_voila_lsv(gene_id, lsv_id) for gene_id, lsv_id in lsvs)
#                     gene_list.append(gene)
#
#                 if len(gene_list) == constants.MAX_GENES:
#                     yield lsv_dict, gene_list
#                     gene_list = []
#                     lsv_dict = {}
#
#             if gene_list:
#                 yield lsv_dict, gene_list
#
#     def get_gene(self, gene_id):
#         """
#         Get gene by its gene id.
#         :param gene_id: unique gene id
#         :return: GeneGraphics
#         """
#         genes = self.hdf5[self.GENES]
#         try:
#             gene = genes[gene_id]
#         except KeyError:
#             raise GeneIdNotFoundInSpliceGraphFile(gene_id)
#
#         return GeneGraphic.easy_from_hdf5(gene)
#
#     def add_experiment_names(self, experiment_names):
#         """
#         Add experiment names to splice graph.
#         :param experiment_names: list of experiment names
#         :return: None
#         """
#         # self.hdf5[self.ROOT].attrs[EXPERIMENT_NAMES] = list(experiment_names)
#         HDF5.create(self.hdf5[self.ROOT].attrs, EXPERIMENT_NAMES, experiment_names)
#
#     def get_experiments(self):
#         """
#         Get list of experiment names from splice graph.
#         :return: list
#         """
#         return self.hdf5[self.ROOT].attrs[EXPERIMENT_NAMES]
#
#     def check_version(self):
#         if self.file_version != constants.SPLICE_GRAPH_FILE_VERSION:
#             voila_log().warning('Splice graph file version isn\'t current.  This will probably cause significant '
#                                 'issues with the voila output.  It would be best to run build again with the current '
#                                 'version of MAJIQ.')


class SpliceGraph(Genes, Junctions, Exons):
    pass


class Voila(VoilaHDF5):
    pass
