import os

from sqlalchemy import create_engine, event, exists, and_
from sqlalchemy.orm import sessionmaker

from voila.api import splice_graph_model as model

Session = sessionmaker()


class SpliceGraphSQL():
    def __init__(self, filename, delete=False):
        self.filename = filename
        if delete is True:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

        engine = create_engine('sqlite:///{0}'.format(filename))
        event.listen(engine, 'connect', self._fk_pragma_on_connect)
        model.Base.metadata.create_all(engine)
        Session.configure(bind=engine)
        self.session = Session()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @staticmethod
    def _fk_pragma_on_connect(dbapi_con, con_record):
        dbapi_con.execute('pragma foreign_keys=ON')

    def close(self):
        self.commit()
        self.session.close_all()

    def add_experiment_names(self, experiment_names):
        session = self.session
        session.add_all([model.Experiment(name=name) for name in experiment_names])

    def get_experiment_names(self):
        return (e for e, in self.session.query(model.Experiment.name).all())

    def commit(self):
        self.session.commit()


class Exons(SpliceGraphSQL):
    class _Exon:
        def __init__(self, session):
            self.session = session

        def add(self, gene_id, start, end, **kwargs):
            coords_extra = kwargs.pop('coords_extra', [])
            alt_ends = kwargs.pop('alt_ends', [])
            alt_starts = kwargs.pop('alt_starts', [])

            exon = model.Exon(gene_id=gene_id, start=start, end=end, **kwargs)

            for ce_start, ce_end in coords_extra:
                exon.coords_extra.append(model.CoordsExtra(start=int(ce_start), end=int(ce_end)))
            for alt_end in alt_ends:
                exon.alt_ends.append(model.AltEnds(coordinate=int(alt_end)))
            for alt_start in alt_starts:
                exon.alt_starts.append(model.AltStarts(coordinate=int(alt_start)))

            self.session.add(exon)

            return exon

        def update(self, gene_id, start, end, **kwargs):
            pass

        def get(self, gene_id, start, end):
            # return self.session.query(model.Exon).get('{0}:{1}-{2}'.format(gene_id, start, end))
            return self.session.query(model.Exon).get((gene_id, start, end))

    @property
    def exon(self):
        return self._Exon(self.session)

    @property
    def exons(self):
        return self.session.query(model.Exon).all()


class Junctions(SpliceGraphSQL):
    class _Junction:
        def __init__(self, session):
            self.session = session

        def add(self, gene_id, start, end, **kwargs):
            reads = kwargs.pop('reads', [])

            junc = model.Junction(gene_id=gene_id, start=start, end=end, **kwargs)

            for r, e in reads:
                junc.reads.append(model.Reads(reads=int(r), experiment_id=e))

            self.session.add(junc)

            return junc

        def exists(self, gene_id, start, end):
            return self.session.query(
                exists().where(and_(model.Junction.gene_id == gene_id, model.Junction.start == start,
                                    model.Junction.end == end))).scalar()

        def update(self, gene_id, start, end, **kwargs):
            reads = kwargs.pop('reads', [])
            j = self.get(gene_id, start, end)
            for k, v in kwargs.items():
                setattr(j, k, v)
            for r, e in reads:
                j.reads.append(model.Reads(reads=int(r), experiment_id=e))

        def get(self, gene_id, start, end):
            return self.session.query(model.Junction).get((gene_id, start, end))

        def update_reads(self, gene_id, start, end, reads, experiment):
            r = model.Reads(junction_gene_id=gene_id, junction_start=start, junction_end=end,
                            experiment_name=experiment, reads=reads)
            self.session.add(r)

    @property
    def junction(self):
        return self._Junction(self.session)

    @property
    def junctions(self):
        return self.session.query(model.Junction).all()


class Genes(SpliceGraphSQL):
    class _Gene:
        def __init__(self, session):
            self.session = session

        def add(self, gene_id, **kwargs):
            g = model.Gene(id=gene_id, **kwargs)
            self.session.add(g)
            return g

        def get(self, gene_id):
            return self.session.query(model.Gene).get(gene_id)

        def update(self):
            pass

    @property
    def genes(self):
        return self.session.query(model.Gene).all()

    @property
    def gene(self):
        return self._Gene(self.session)
