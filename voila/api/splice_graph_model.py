from sqlalchemy import Integer, String, Column, Boolean, ForeignKey, ForeignKeyConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from voila import constants

Base = declarative_base()


def base_iter(self):
    for col in self.__table__.columns.keys():
        yield col, getattr(self, col)


class Experiment(Base):
    __tablename__ = 'experiment'
    __iter__ = base_iter
    name = Column(String, primary_key=True)


class CoordsExtra(Base):
    __tablename__ = 'coords_extra'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    start = Column(Integer)
    end = Column(Integer)
    exon_gene_id = Column(String)
    exon_start = Column(Integer)
    exon_end = Column(Integer)

    __table_args__ = (
        ForeignKeyConstraint([exon_gene_id, exon_start, exon_end], ['exon.gene_id', 'exon.start', 'exon.end']),)


class AltStarts(Base):
    __tablename__ = 'alt_starts'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    coordinate = Column(Integer)
    exon_gene_id = Column(String)
    exon_start = Column(Integer)
    exon_end = Column(Integer)

    __table_args__ = (
        ForeignKeyConstraint([exon_gene_id, exon_start, exon_end], ['exon.gene_id', 'exon.start', 'exon.end']),)


class AltEnds(Base):
    __tablename__ = 'alt_ends'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    coordinate = Column(Integer)
    exon_gene_id = Column(String)
    exon_start = Column(Integer)
    exon_end = Column(Integer)

    __table_args__ = (
        ForeignKeyConstraint([exon_gene_id, exon_start, exon_end], ['exon.gene_id', 'exon.start', 'exon.end']),)


class Reads(Base):
    __tablename__ = 'reads'
    __iter__ = base_iter

    reads = Column(Integer, nullable=False)

    experiment_name = Column(String, ForeignKey('experiment.name'), primary_key=True)
    junction_gene_id = Column(String, primary_key=True)
    junction_start = Column(Integer, primary_key=True)
    junction_end = Column(Integer, primary_key=True)

    __table_args__ = (
        ForeignKeyConstraint([junction_gene_id, junction_start, junction_end],
                             ['junction.gene_id', 'junction.start', 'junction.end']),)

    experiment = relationship('Experiment')
    junction = relationship('Junction')


class Junction(Base):
    __tablename__ = 'junction'
    __iter__ = base_iter
    gene_id = Column(String, ForeignKey('gene.id'), primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)

    intron_retention = Column(Integer)
    annotated = Column(Boolean)

    reads = relationship('Reads')
    gene = relationship('Gene')

    @property
    def length(self):
        return self.end - self.start

    def get_junction_type(self, experiment_name):
        reads = self.get_reads(experiment_name)

        if self.annotated and reads == 0:
            reads_sum = sum(r.reads for r in self.reads)
            if reads_sum - reads > 0:
                return constants.JUNCTION_TYPE_DB_OTHER_RNASEQ
            else:
                return constants.JUNCTION_TYPE_DB

        if self.annotated and reads > 0:
            return constants.JUNCTION_TYPE_DB_RNASEQ

        if not self.annotated and reads > 0:
            return constants.JUNCTION_TYPE_RNASEQ

        return constants.JUNCTION_TYPE_RNASEQ

    def get_reads(self, experiment_name):
        return next((r.reads for r in self.reads if r.experiment_name == experiment_name), 0)

    def get_experiment(self, experiment_name):
        junction = dict(self)

        del junction['gene_id']

        junction['reads'] = self.get_reads(experiment_name)
        junction['junction_type'] = self.get_junction_type(experiment_name)
        junction['length'] = self.length

        if self.intron_retention:
            junction['intron_retention'] = -1
            for exon in (e for e in self.gene.exons if not e.intron_retention):
                if self.start == exon.end:
                    junction['intron_retention'] = constants.IR_TYPE_START
                elif self.end == exon.start:
                    junction['intron_retention'] = constants.IR_TYPE_END

        return junction


class Exon(Base):
    __tablename__ = 'exon'
    __iter__ = base_iter

    gene_id = Column(String, ForeignKey('gene.id'), primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)

    intron_retention = Column(Integer)
    annotated = Column(Boolean)

    coords_extra = relationship('CoordsExtra')
    alt_starts = relationship('AltStarts')
    alt_ends = relationship('AltEnds')
    gene = relationship('Gene')

    @property
    def a5(self):
        def a5_filter(js):
            for j in js:
                if self.end == -1:
                    continue

                if self.start == -1 and j.start == self.end:
                    yield j

                if self.start <= j.start <= self.end:
                    yield j

        junctions = sorted(self.gene.junctions, key=lambda j: (j.start, j.end))
        return tuple(junctions.index(j) for j in a5_filter(junctions))

    @property
    def a3(self):
        def a3_filter(js):
            for j in js:
                if self.start == -1:
                    continue

                if self.end == -1 and j.end == self.start:
                    yield j

                if self.start <= j.end <= self.end:
                    yield j

        junctions = sorted(self.gene.junctions, key=lambda j: (j.start, j.end))
        return tuple(junctions.index(j) for j in a3_filter(junctions))

    @property
    def length(self):
        if -1 in (self.start, self.end):
            return 10
        return self.end - self.start

    def get_exon_type(self, experiment_name):

        if self.start == -1:
            return constants.EXON_TYPE_MISSING_START
        if self.end == -1:
            return constants.EXON_TYPE_MISSING_END

        junctions = sorted(self.gene.junctions, key=lambda j: (j.start, j.end))
        has_reads = any(junctions[x].get_reads(experiment_name) for x in self.a5 + self.a3)

        if self.annotated and not has_reads:
            return constants.EXON_TYPE_DB

        if self.annotated and has_reads:
            return constants.EXON_TYPE_DB_RNASEQ

        if not self.annotated and has_reads:
            return constants.EXON_TYPE_RNASEQ

        return constants.EXON_TYPE_RNASEQ

    def get_experiment(self, experiment_name):
        exon = dict(self)
        exon['a5'] = self.a5
        exon['a3'] = self.a3
        exon['length'] = self.length
        exon['exon_type'] = self.get_exon_type(experiment_name)

        if self.start == -1:
            exon['start'] = self.end - 10
        elif self.end == -1:
            exon['end'] = self.start + 10

        del exon['gene_id']

        return exon


class Gene(Base):
    __tablename__ = 'gene'
    __iter__ = base_iter

    id = Column(String, primary_key=True)
    name = Column(String)
    strand = Column(String)
    chromosome = Column(String)

    junctions = relationship('Junction')
    exons = relationship('Exon')

    @property
    def start(self):
        return sorted(e.start for e in self.exons if e.start != -1)[0]

    @property
    def end(self):
        return sorted((e.end for e in self.exons if e.end != -1), reverse=True)[0]

    def get_experiment(self, experiment_name):
        gene = dict(self)
        gene['_id'] = self.id + '_' + experiment_name

        exons = tuple(e.get_experiment(experiment_name) for e in self.exons)
        exons = sorted(exons, key=lambda e: (e['start'], e['end']))
        gene['exons'] = exons

        juncs = sorted(self.junctions, key=lambda j: (j.start, j.end))
        gene['junctions'] = tuple(j.get_experiment(experiment_name) for j in juncs)

        gene['start'] = gene['exons'][0]['start']
        gene['end'] = gene['exons'][-1]['end']
        gene['length'] = gene['end'] - gene['start']

        return gene
