from sqlalchemy import Integer, String, Column, Boolean, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


def base_iter(self):
    for col in self.__table__.columns.keys():
        yield col, getattr(self, col)


class Experiment(Base):
    __tablename__ = 'experiment'
    __iter__ = base_iter
    id = Column(String, primary_key=True)


class CoordsExtra(Base):
    __tablename__ = 'coords_extra'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    start = Column(Integer)
    end = Column(Integer)
    exon_id = Column(String, ForeignKey('exon.id'))


class AltStarts(Base):
    __tablename__ = 'alt_starts'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    coordinate = Column(Integer)
    exon_id = Column(String, ForeignKey('exon.id'))


class AltEnds(Base):
    __tablename__ = 'alt_ends'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    coordinate = Column(Integer)
    exon_id = Column(String, ForeignKey('exon.id'))


class Reads(Base):
    __tablename__ = 'reads'
    __iter__ = base_iter
    id = Column(Integer, primary_key=True)
    reads = Column(Integer)
    experiment_id = Column(String, ForeignKey('experiment.id'))
    experiment = relationship('Experiment')
    junction_id = Column(String, ForeignKey('junction.id'))
    junction = relationship('Junction')


class Junction(Base):
    __tablename__ = 'junction'
    __iter__ = base_iter
    id = Column(String, primary_key=True)
    reads = relationship('Reads')
    intron_retention = Column(Integer)
    annotated = Column(Boolean)
    gene_id = Column(String, ForeignKey('gene.id'))
    gene = relationship('Gene')

    @property
    def start(self):
        return self.id.split(':')[1].split('-')[0]

    @property
    def end(self):
        return self.id.split(':')[1].split('-')[1]


class Exon(Base):
    __tablename__ = 'exon'
    __iter__ = base_iter
    id = Column(String, primary_key=True)
    coords_extra = relationship('CoordsExtra')
    intron_retention = Column(Integer)
    alt_starts = relationship('AltStarts')
    alt_ends = relationship('AltEnds')
    annotated = Column(Boolean)
    gene_id = Column(String, ForeignKey('gene.id'))

    @property
    def start(self):
        return self.id.split(':')[1].split('-')[0]

    @property
    def end(self):
        return self.id.split(':')[1].split('-')[1]


class Gene(Base):
    __tablename__ = 'gene'
    __iter__ = base_iter
    id = Column(String, primary_key=True)
    name = Column(String)
    strand = Column(String)
    chromosome = Column(String)
    junctions = relationship('Junction')
    exons = relationship('Exon')
