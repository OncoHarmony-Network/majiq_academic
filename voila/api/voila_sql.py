from sqlalchemy import exists, and_

from voila.api import voila_model as model
from voila.api.sql import SQL, SQLType


class VoilaSQL(SQL):
    def __init__(self, filename, delete=False):
        super().__init__(filename, delete=delete, model=model)
        self.metadata = model.Metadata()
        self.session.add(self.metadata)

    def add_genome(self, genome):
        self.metadata.genome = genome

    def set_analysis_type(self, analysis_type):
        self.metadata.analysis_type = analysis_type

    def add_experiments(self, group_name, experiment_names):
        self.metadata.experiment_names = experiment_names
        self.metadata.group_name = group_name


class Exons(VoilaSQL):
    class _Exons(SQLType):
        def __init__(self, session, sql, lsv_id, start, end):
            super().__init__(session)
            self.lsv_id = lsv_id
            self.start = start
            self.end = end
            self.sql = sql

        @property
        def exists(self):
            return self.session.query(
                exists().where(and_(model.Exon.lsv_id == self.lsv_id, model.Exon.start == self.start,
                                    model.Exon.end == self.end))).scalar()

        def add(self, **kwargs):
            e = model.Exon(lsv_id=self.lsv_id, start=self.start, end=self.end, **kwargs)
            self.session.add(e)
            return e

        @property
        def get(self):
            return self.session.query(model.Lsv).get((self.lsv_id, self.start, self.end))

    @property
    def exons(self):
        return self.session.query(model.Exon).all()

    def exon(self, lsv_id, start, end):
        return self._Exons(self.session, self, lsv_id, start, end)


class Junctions(VoilaSQL):
    class _Junctions(SQLType):
        def __init__(self, session, lsv_id, start, end):
            super().__init__(session)
            self.lsv_id = lsv_id
            self.start = start
            self.end = end

        @property
        def exists(self):
            return self.session.query(
                exists().where(and_(model.Junction.lsv_id == self.lsv_id, model.Junction.start == self.start,
                                    model.Junction.end == self.end))).scalar()

        def add(self, **kwargs):
            j = model.Junction(lsv_id=self.lsv_id, start=self.start, end=self.end, **kwargs)
            self.session.add(j)
            return j

        @property
        def get(self):
            return self.session.query(model.Junction).get((self.lsv_id, self.start, self.end))

    @property
    def junctions(self):
        return self.session.query(model.Junction).all()

    def junction(self, lsv_id, start, end):
        return self._Junctions(self.session, lsv_id, start, end)


class Lsvs(VoilaSQL):
    class _(SQLType):
        def __init__(self, session, lsv_id):
            super().__init__(session)
            self.lsv_id = lsv_id

        def add(self, **kwargs):
            l = model.Lsv(id=self.lsv_id, **kwargs)
            self.session.add(l)
            return l

        @property
        def get(self, *args):
            return self.session.query(model.Lsv).get(self.lsv_id)

        @property
        def exists(self):
            return self.session.query(exists().where(model.Lsv.id == self.lsv_id)).scalar()

    @property
    def lsvs(self):
        return self.session.query(model.Lsv).all()

    def lsv(self, lsv_id):
        return self._(self.session, lsv_id)


class Voila(Lsvs, Junctions, Exons):
    pass
