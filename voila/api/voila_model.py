import pickle

from sqlalchemy import Column, String, Integer, ForeignKey, Binary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


class Metadata(Base):
    __tablename__ = 'metadata'
    id = Column(Integer, primary_key=True)
    genome = Column(String)


class Lsv(Base):
    __tablename__ = 'lsv'
    id = Column(String, primary_key=True)
    type = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    gene_name = Column(String)
    strand = Column(String)
    chromosome = Column(String)

    exons = relationship('Exon')
    junctions = relationship('Junction')


class Exon(Base):
    __tablename__ = 'exon'
    lsv_id = Column(String, ForeignKey('lsv.id'), primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)

    lsv = relationship('Lsv')


class pickled_property(property):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.attr_name = '_{0}'.format(args[0].__name__)

    def __get__(self, *args, **kwargs):
        obj = args[0]
        if obj is not None:
            return pickle.loads(getattr(obj, self.attr_name))

    def __set__(self, *args, **kwargs):
        obj, value = args
        setattr(obj, self.attr_name, pickle.dumps(list(map(float, value))))


class Junction(Base):
    __tablename__ = 'junction'
    lsv_id = Column(String, ForeignKey('lsv.id'), primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)

    _bin = Column(Binary)
    _means = Column(Binary)
    _psi1 = Column(Binary)
    _psi2 = Column(Binary)
    _means_psi1 = Column(Binary)
    _means_psi2 = Column(Binary)

    lsv = relationship('Lsv')

    @pickled_property
    def bin(self):
        pass

    @pickled_property
    def means(self):
        pass

    @pickled_property
    def psi1(self):
        pass

    @pickled_property
    def psi2(self):
        pass

    @pickled_property
    def means_psi1(self):
        pass

    @pickled_property
    def means_psi2(self):
        pass
