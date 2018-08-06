from abc import ABC, abstractmethod


class SpliceGraphSQLAbstract(ABC):

    @property
    @abstractmethod
    def genome(self):
        pass

    @genome.setter
    @abstractmethod
    def genome(self, g):
        pass

    @property
    @abstractmethod
    def experiment_names(self):
        pass

    @experiment_names.setter
    @abstractmethod
    def experiment_names(self, names):
        pass

    @property
    @abstractmethod
    def file_version(self):
        pass

    @file_version.setter
    @abstractmethod
    def file_version(self, version):
        pass


class SpliceGraphType(ABC):
    def __bool__(self):
        return self.exists

    def __iter__(self):
        return self.get.__iter__()

    @abstractmethod
    def add(self, *args, **kwargs):
        pass

    @property
    @abstractmethod
    def get(self):
        pass

    @property
    @abstractmethod
    def exists(self):
        pass


class GenesAbstract(ABC):
    @property
    @abstractmethod
    def genes(self):
        pass

    @abstractmethod
    def gene(self, gene_id):
        pass


class ExonsAbstract(ABC):
    @abstractmethod
    def exon(self, gene_id: str, start: int, end: int):
        pass

    @property
    @abstractmethod
    def exons(self):
        pass


class JunctionsAbstract(ABC):
    @abstractmethod
    def junction(self, gene_id: str, start: int, end: int):
        pass

    @property
    @abstractmethod
    def junctions(self):
        pass


class IntronRetentionAbstract(ABC):
    @abstractmethod
    def intron_retention(self, gene_id: str, start: int, end: int):
        pass

    @property
    @abstractmethod
    def intron_retentions(self):
        pass
