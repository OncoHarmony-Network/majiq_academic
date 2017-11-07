from sqlalchemy import exists, and_
from sqlalchemy.orm import sessionmaker

from voila.api import splice_graph_model as model
from voila.api.sql import SQL, SQLType

Session = sessionmaker()


class SpliceGraphSQL(SQL):
    def __init__(self, filename, delete=False):
        super().__init__(filename, delete=delete, model=model)

    def add_experiment_names(self, experiment_names):
        session = self.session
        session.add_all([model.Experiment(name=name) for name in experiment_names])

    def get_experiment_names(self):
        return (e for e, in self.session.query(model.Experiment.name).all())


class Exons(SpliceGraphSQL):
    class _Exon(SQLType):
        def __init__(self, session, gene_id, start, end):
            super().__init__(session)
            self.gene_id = gene_id
            self.start = start
            self.end = end

        def add(self, **kwargs):
            coords_extra = kwargs.pop('coords_extra', [])
            alt_ends = kwargs.pop('alt_ends', [])
            alt_starts = kwargs.pop('alt_starts', [])

            exon = model.Exon(gene_id=self.gene_id, start=self.start, end=self.end, **kwargs)

            for ce_start, ce_end in coords_extra:
                exon.coords_extra.append(model.CoordsExtra(start=int(ce_start), end=int(ce_end)))
            for alt_end in alt_ends:
                exon.alt_ends.append(model.AltEnds(coordinate=int(alt_end)))
            for alt_start in alt_starts:
                exon.alt_starts.append(model.AltStarts(coordinate=int(alt_start)))

            self.session.add(exon)

            return exon

        @property
        def get(self):
            return self.session.query(model.Exon).get((self.gene_id, self.start, self.end))

        @property
        def exists(self):
            return self.session.query(
                exists().where(and_(model.Exon.gene_id == self.gene_id, model.Exon.start == int(self.start),
                                    model.Exon.end == int(self.end)))).scalar()

    def exon(self, gene_id, start, end):
        return self._Exon(self.session, gene_id, start, end)

    @property
    def exons(self):
        return self.session.query(model.Exon).all()


class Junctions(SpliceGraphSQL):
    class _Junction(SQLType):
        def __init__(self, session, gene_id, start, end):
            super().__init__(session)
            self.gene_id = gene_id
            self.start = int(start)
            self.end = int(end)

        def add(self, **kwargs):
            reads = kwargs.pop('reads', [])

            junc = model.Junction(gene_id=self.gene_id, start=self.start, end=self.end, **kwargs)

            for r, e in reads:
                junc.reads.append(model.Reads(reads=int(r), experiment_id=e))

            self.session.add(junc)

            return junc

        def exists(self):
            return self.session.query(
                exists().where(and_(model.Junction.gene_id == self.gene_id, model.Junction.start == self.start,
                                    model.Junction.end == self.end))).scalar()

        def get(self):
            return self.session.query(model.Junction).get((self.gene_id, self.start, self.end))

        def update_reads(self, reads, experiment):
            r = model.Reads(junction_gene_id=self.gene_id, junction_start=self.start, junction_end=self.end,
                            experiment_name=experiment, reads=int(reads))
            self.session.add(r)

    def junction(self, gene_id, start, end):
        return self._Junction(self.session, gene_id, start, end)

    @property
    def junctions(self):
        return self.session.query(model.Junction).all()


class Genes(SpliceGraphSQL):
    class _Gene(SQLType):
        def __init__(self, session, gene_id):
            super().__init__(session)
            self.gene_id = gene_id

        def add(self, **kwargs):
            g = model.Gene(id=self.gene_id, **kwargs)
            self.session.add(g)
            return g

        def get(self):
            return self.session.query(model.Gene).get(self.gene_id)

        def exists(self):
            return self.session.query(exists().where(model.Gene.id == self.gene_id)).scalar()

    @property
    def genes(self):
        return self.session.query(model.Gene).all()

    def gene(self, gene_id):
        return self._Gene(self.session, gene_id)
