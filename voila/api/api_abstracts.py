from abc import ABC, abstractmethod


class VoilaAbstract(ABC):
    @abstractmethod
    def __enter__(self):
        pass

    @abstractmethod
    def __exit__(self, type, value, traceback):
        pass

    @abstractmethod
    def _metainfo(self):
        pass

    @abstractmethod
    def close(self):
        pass

    @abstractmethod
    def check_version(self):
        pass

    @abstractmethod
    def get_metainfo(self):
        pass

    @abstractmethod
    def get_voila_lsv(self, gene_id, lsv_id):
        pass

    @abstractmethod
    def get_lsv_count(self, args):
        pass

    @abstractmethod
    def get_voila_lsvs(self, args=None, gene_id=None):
        pass

    @abstractmethod
    def get_gene_ids(self, args=None):
        pass

    @abstractmethod
    def get_lsvs(self, args=None, gene_id=None):
        pass

    @abstractmethod
    def get_lsv_means(self, gene_id, lsv_id):
        pass

    @abstractmethod
    def get_lsv(self, gene_id, lsv_id):
        pass

    @abstractmethod
    def get_lsv_type(self, gene_id, lsv_id):
        pass

    @abstractmethod
    def get_gene_name(self, gene_id, lsv_id):
        pass

    @abstractmethod
    def get_lsv_ids(self, gene_id):
        pass

    @abstractmethod
    def add_stat_names(self, stat_names):
        pass

    @abstractmethod
    def add_lsv(self, voilaLsv):
        pass

    @abstractmethod
    def add_genome(self, genome):
        pass

    @abstractmethod
    def add_metainfo(self, genome, group1, experiments1, group2=None, experiments2=None):
        pass

    @abstractmethod
    def add_experiments(self, group_name, experiment_names):
        pass

    @abstractmethod
    def set_analysis_type(self, analysis_type):
        pass


class SpliceGraphsAbstract(ABC):
    @abstractmethod
    def __enter__(self):
        pass

    @abstractmethod
    def __exit__(self, type, value, traceback):
        pass

    @abstractmethod
    def close(self):
        pass

    @abstractmethod
    def erase_splice_graph_file(self):
        pass

    @abstractmethod
    def check_version(self):
        pass

    @abstractmethod
    def get_page_count(self, args):
        pass

    @abstractmethod
    def get_gene_name(self, gene_id):
        pass

    @abstractmethod
    def get_gene_ids(self, args=None):
        pass

    @abstractmethod
    def get_genes(self):
        pass

    @abstractmethod
    def get_paginated_genes(self, args):
        pass

    @abstractmethod
    def get_paginated_genes_with_lsvs(self, args):
        pass

    @abstractmethod
    def get_gene(self, gene_id):
        pass

    @abstractmethod
    def get_experiments(self):
        pass

    @abstractmethod
    def add_gene(self, gene):
        pass

    @abstractmethod
    def add_experiment_names(self, experiment_names):
        pass
