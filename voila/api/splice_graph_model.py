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


class AltStart(Base):
    __tablename__ = 'alt_start'
    __iter__ = base_iter
    gene_id = Column(String, primary_key=True)
    coordinate = Column(Integer, primary_key=True)

    __table_args__ = (ForeignKeyConstraint((gene_id,), ['gene.id']),)


class AltEnd(Base):
    __tablename__ = 'alt_end'
    __iter__ = base_iter
    gene_id = Column(String, primary_key=True)
    coordinate = Column(Integer, primary_key=True)

    __table_args__ = (ForeignKeyConstraint((gene_id,), ['gene.id']),)


class JunctionReads(Base):
    __tablename__ = 'junction_reads'
    __iter__ = base_iter

    reads = Column(Integer, nullable=False)
    experiment_name = Column(String, ForeignKey('experiment.name'), primary_key=True)
    junction_gene_id = Column(String, primary_key=True)
    junction_start = Column(Integer, primary_key=True)
    junction_end = Column(Integer, primary_key=True)

    __table_args__ = (
        ForeignKeyConstraint((junction_gene_id, junction_start, junction_end),
                             ['junction.gene_id', 'junction.start', 'junction.end']),)


class IntronRetentionReads(Base):
    __tablename__ = 'intron_retention_reads'
    __iter__ = base_iter

    reads = Column(Integer, nullable=False)
    experiment_name = Column(String, ForeignKey('experiment.name'), primary_key=True)
    intron_retention_gene_id = Column(String, primary_key=True)
    intron_retention_start = Column(Integer, primary_key=True)
    intron_retention_end = Column(Integer, primary_key=True)

    __table_args__ = (
        ForeignKeyConstraint((intron_retention_gene_id, intron_retention_start, intron_retention_end),
                             ['intron_retention.gene_id', 'intron_retention.start', 'intron_retention.end']),)


class Junction(Base):
    __tablename__ = 'junction'
    __iter__ = base_iter
    gene_id = Column(String, ForeignKey('gene.id'), primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)
    has_reads = Column(Boolean, default=False)
    annotated = Column(Boolean)

    gene = relationship('Gene')
    reads = relationship('JunctionReads')


class Exon(Base):
    __tablename__ = 'exon'
    __iter__ = base_iter
    __contains__ = exon_contains

    gene_id = Column(String, ForeignKey('gene.id'), primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)
    annotated_start = Column(Integer, default=-1)
    annotated_end = Column(Integer, default=-1)
    annotated = Column(Boolean)

    gene = relationship('Gene')

    @property
    def a3(self):
        def a3_filter(js):
            if self.end == -1:
                return

            for j in js:
                if j.start in self:
                    yield j

        yield from a3_filter(self.gene.junctions)

    @property
    def a5(self):
        def a5_filter(js):
            if self.start == -1:
                return

            for j in js:
                if j.end in self:
                    yield j

        yield from a5_filter(self.gene.junctions)


class IntronRetention(Base):
    __tablename__ = 'intron_retention'
    __iter__ = base_iter

    gene_id = Column(String, ForeignKey('gene.id'), primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)
    has_reads = Column(Boolean, default=False)
    annotated = Column(Boolean)

    reads = relationship('IntronRetentionReads')
    gene = relationship('Gene')


class Gene(Base):
    __tablename__ = 'gene'
    __iter__ = base_iter

    id = Column(String, primary_key=True)
    name = Column(String)
    strand = Column(String)
    chromosome = Column(String)

    junctions = relationship('Junction')
    intron_retentions = relationship('IntronRetention')
    alt_starts = relationship('AltStart')
    alt_ends = relationship('AltEnd')

    @property
    def start(self):
        return sorted(e.start for e in self.exons if e.start != -1)[0]

    @property
    def end(self):
        return sorted((e.end for e in self.exons if e.end != -1), reverse=True)[0]
