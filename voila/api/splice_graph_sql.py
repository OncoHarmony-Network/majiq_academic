from abc import ABC, abstractmethod
from typing import List

from sqlalchemy import exists, and_
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound

from voila.api import splice_graph_model as model
from voila.api.sql import SQL

Session = sessionmaker()

default_commit_on_count = 1e4


class SpliceGraphSQL(SQL):
    def __init__(self, filename: str, delete: bool = False, nprocs: int = 1):
        super().__init__(filename, model, delete, nprocs)

    @property
    def genome(self):
        return self.session.query(model.Genome.name).one()[0]

    @genome.setter
    def genome(self, g: str):
        self.session.add(model.Genome(name=g))

    @property
    def experiment_names(self):
        return (e for e, in self.session.query(model.Experiment.name).all())

    @experiment_names.setter
    def experiment_names(self, names: List[str]):
        self.session.add_all([model.Experiment(name=name) for name in names])

    @property
    def file_version(self):
        try:
            return self.session.query(model.FileVersion.value).one()[0]
        # when file doesn't contain FileVersion table
        except OperationalError:
            return -1
        # no value in table
        except NoResultFound:
            return -1

    @file_version.setter
    def file_version(self, version: int):
        self.session.add(model.FileVersion(value=version))


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


class IntronRetention(SpliceGraphSQL):
    class _IntronRetention(SpliceGraphType):
        def __init__(self, sql, gene_id: str, start: int, end: int):
            self.sql = sql
            self.gene_id = str(gene_id)
            self.start = int(start)
            self.end = int(end)

        def add(self, **kwargs):
            ir = model.IntronRetention(gene_id=self.gene_id, start=self.start, end=self.end, **kwargs)
            self.sql.add(ir)
            self.sql.commit(default_commit_on_count)
            return ir

        @property
        def get(self):
            return self.sql.session.query(model.IntronRetention).get((self.gene_id, self.start, self.end))

        @property
        def exists(self):
            return self.sql.session.query(
                exists().where(
                    and_(model.IntronRetention.gene_id == self.gene_id, model.IntronRetention.start == self.start,
                         model.IntronRetention.end == self.end))).scalar()

        def reads_exists(self, experiment, reads):
            experiment = str(experiment)
            reads = int(reads)
            return self.sql.session.query(
                exists().where(and_(
                    model.IntronRetentionReads.intron_retention_gene_id == self.gene_id,
                    model.IntronRetentionReads.intron_retention_start == self.start,
                    model.IntronRetentionReads.intron_retention_end == self.end,
                    model.IntronRetentionReads.experiment_name == experiment,
                    model.IntronRetentionReads.reads == reads,
                ))).scalar()

        def update_reads(self, experiment: str, reads: int):
            experiment = str(experiment)
            reads = int(reads)

            if reads:
                r = model.IntronRetentionReads(intron_retention_gene_id=self.gene_id, intron_retention_start=self.start,
                                               intron_retention_end=self.end, experiment_name=experiment, reads=reads)

                with self.sql.session.no_autoflush:
                    if self.reads_exists(experiment, reads):
                        return -1
                    else:
                        self.get.has_reads = True
                        self.sql.add(r)
                        self.sql.commit(default_commit_on_count)

            return 0

    def intron_retention(self, gene_id: str, start: int, end: int):
        return self._IntronRetention(self, gene_id, start, end)

    @property
    def intron_retentions(self):
        return self.session.query(model.IntronRetention).all()


class Exons(SpliceGraphSQL):
    class _Exon(SpliceGraphType):
        def __init__(self, sql, gene_id: str, start: int, end: int):
            self.sql = sql
            self.gene_id = str(gene_id)
            self.start = int(start)
            self.end = int(end)

        def add(self, **kwargs):
            coords_extra = kwargs.pop('coords_extra', [])
            alt_ends = kwargs.pop('alt_ends', [])
            alt_starts = kwargs.pop('alt_starts', [])

            exon = model.Exon(gene_id=self.gene_id, start=self.start, end=self.end, **kwargs)
            exon.coords_extra = [model.CoordsExtra(start=int(s), end=int(e)) for s, e in coords_extra]
            exon.alt_ends = [model.AltEnds(coordinate=int(alt_end)) for alt_end in alt_ends]
            exon.alt_starts = [model.AltStarts(coordinate=int(alt_start)) for alt_start in alt_starts]

            self.sql.add(exon)
            self.sql.commit(default_commit_on_count)
            return exon

        @property
        def get(self):
            return self.sql.session.query(model.Exon).get((self.gene_id, self.start, self.end))

        @property
        def exists(self):
            return self.sql.session.query(
                exists().where(and_(model.Exon.gene_id == self.gene_id, model.Exon.start == self.start,
                                    model.Exon.end == self.end))).scalar()

    def exon(self, gene_id: str, start: int, end: int):
        return self._Exon(self, gene_id, start, end)

    @property
    def exons(self):
        return self.session.query(model.Exon).all()


class Junctions(SpliceGraphSQL):
    class _Junction(SpliceGraphType):
        def __init__(self, sql, gene_id: str, start: int, end: int):
            self.sql = sql
            self.gene_id = str(gene_id)
            self.start = int(start)
            self.end = int(end)

        def add(self, **kwargs):
            reads = kwargs.pop('reads', [])

            junc = model.Junction(gene_id=self.gene_id, start=self.start, end=self.end, **kwargs)
            junc.reads = [model.JunctionReads(reads=int(r), experiment_name=e) for r, e in reads]

            self.sql.add(junc)
            self.sql.commit(default_commit_on_count)
            return junc

        @property
        def exists(self):
            return self.sql.session.query(
                exists().where(and_(model.Junction.gene_id == self.gene_id, model.Junction.start == self.start,
                                    model.Junction.end == self.end))).scalar()

        @property
        def get(self):
            return self.sql.session.query(model.Junction).get((self.gene_id, self.start, self.end))

        def update_reads(self, experiment: str, reads: int):
            reads = int(reads)
            experiment = str(experiment)

            if not reads:
                return

            r = model.JunctionReads(junction_gene_id=self.gene_id, junction_start=self.start, junction_end=self.end,
                                    experiment_name=experiment, reads=reads)

            with self.sql.session.no_autoflush:
                self.get.has_reads = True
                self.sql.add(r)
                self.sql.commit(default_commit_on_count)

    def junction(self, gene_id: str, start: int, end: int):
        return self._Junction(self, gene_id, start, end)

    @property
    def junctions(self):
        return self.session.query(model.Junction).all()


class Genes(SpliceGraphSQL):
    class _Gene(SpliceGraphType):
        def __init__(self, sql, gene_id: str):
            self.sql = sql
            self.gene_id = str(gene_id)

        def add(self, **kwargs):
            g = model.Gene(id=self.gene_id, **kwargs)
            self.sql.add(g)
            self.sql.commit(default_commit_on_count)
            return g

        @property
        def get(self):
            return self.sql.session.query(model.Gene).get(self.gene_id)

        @property
        def exists(self):
            return self.sql.session.query(exists().where(model.Gene.id == self.gene_id)).scalar()

    @property
    def genes(self):
        return self.session.query(model.Gene).all()

    def gene(self, gene_id: str):
        return self._Gene(self, gene_id)
