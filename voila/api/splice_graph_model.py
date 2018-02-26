from sqlalchemy import Integer, String, Column, Boolean, ForeignKey, ForeignKeyConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


def exon_contains(self, coord):
    if -1 in (self.start, self.end):
        return coord in (self.start, self.end)
    else:
        return self.start <= coord <= self.end


def base_iter(self):
    for col in self.__table__.columns.keys():
        yield col, getattr(self, col)


class Genome(Base):
    __tablename__ = 'genome'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    name = Column(String)


class FileVersion(Base):
    __tablename__ = 'file_version'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    value = Column(Integer)


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


class ReadsSum(Base):
    __tablename__ = 'reads_sum'
    __iter__ = base_iter

    sum = Column(Integer, nullable=False)
    junction_gene_id = Column(String, primary_key=True)
    junction_start = Column(Integer, primary_key=True)
    junction_end = Column(Integer, primary_key=True)

    __table_args__ = (
        ForeignKeyConstraint([junction_gene_id, junction_start, junction_end],
                             ['junction.gene_id', 'junction.start', 'junction.end']),)

    junction = relationship('Junction')


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


class Exon(Base):
    __tablename__ = 'exon'
    __iter__ = base_iter
    __contains__ = exon_contains

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
    def a3(self):
        def a3_filter(js):
            if self.end == -1:
                return

            for j in js:
                if j.start in self and j.end not in self:
                    yield j

        yield from a3_filter(self.gene.junctions)

    @property
    def a5(self):
        def a5_filter(js):
            if self.start == -1:
                return

            for j in js:
                if j.end in self and j.start not in self:
                    yield j

        yield from a5_filter(self.gene.junctions)


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
