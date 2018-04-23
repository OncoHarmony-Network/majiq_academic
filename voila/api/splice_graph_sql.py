from abc import ABC, abstractmethod

from sqlalchemy import exists, and_
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound

from voila.api import splice_graph_model as model
from voila.api.sql import SQL

Session = sessionmaker()

default_commit_on_count = 100000


class SpliceGraphSQL(SQL):
    def __init__(self, filename, delete=False):
        super().__init__(filename, model, delete)

    @property
    def genome(self):
        return self.session.query(model.Genome.name).one()[0]

    @genome.setter
    def genome(self, g):
        self.session.add(model.Genome(name=g))

    @property
    def experiment_names(self):
        return (e for e, in self.session.query(model.Experiment.name).all())

    @experiment_names.setter
    def experiment_names(self, names):
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
    def file_version(self, version):
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


class Exons(SpliceGraphSQL):
    class _Exon(SpliceGraphType):
        def __init__(self, sql, gene_id, start, end):
            self.sql = sql
            self.gene_id = gene_id
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
                exists().where(and_(model.Exon.gene_id == self.gene_id, model.Exon.start == int(self.start),
                                    model.Exon.end == int(self.end)))).scalar()

    def exon(self, gene_id, start, end):
        return self._Exon(self, gene_id, start, end)

    @property
    def exons(self):
        return self.session.query(model.Exon).all()


class Junctions(SpliceGraphSQL):
    class _Junction(SpliceGraphType):
        def __init__(self, sql, gene_id, start, end):
            self.sql = sql
            self.gene_id = gene_id
            self.start = int(start)
            self.end = int(end)

        def add(self, **kwargs):
            reads = kwargs.pop('reads', [])

            junc = model.Junction(gene_id=self.gene_id, start=self.start, end=self.end, **kwargs)
            junc.reads = [model.Reads(reads=int(r), experiment_name=e) for r, e in reads]

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

        def update_reads(self, experiment, reads):
            if not int(reads):
                return

            r = model.Reads(junction_gene_id=self.gene_id, junction_start=self.start, junction_end=self.end,
                            experiment_name=experiment, reads=int(reads))
            self.sql.add(r)
            self.get.has_reads = True
            self.sql.commit(default_commit_on_count)

    def junction(self, gene_id, start, end):
        return self._Junction(self, gene_id, start, end)

    @property
    def junctions(self):
        return self.session.query(model.Junction).all()


class Genes(SpliceGraphSQL):
    class _Gene(SpliceGraphType):
        def __init__(self, sql, gene_id):
            self.sql = sql
            self.gene_id = gene_id

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

    def gene(self, gene_id):
        return self._Gene(self, gene_id)
