from sqlalchemy import Integer, String, Column, Boolean, ForeignKey, ForeignKeyConstraint, inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from voila import constants

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

    def get_junction_types(self, experiment_names):
        for experiment_name in experiment_names:
            reads = self.get_reads(experiment_name)

            junc_type = constants.JUNCTION_TYPE_RNASEQ

            if reads == 0 and self.annotated:
                if self.reads_sum(experiment_names) - reads > 0:
                    junc_type = constants.JUNCTION_TYPE_DB_OTHER_RNASEQ
                else:
                    junc_type = constants.JUNCTION_TYPE_DB

            if self.annotated and reads > 0:
                junc_type = constants.JUNCTION_TYPE_DB_RNASEQ

            if not self.annotated and reads > 0:
                junc_type = constants.JUNCTION_TYPE_RNASEQ

            yield experiment_name, junc_type

    def reads_sum(self, experiment_names):
        return sum(self.get_reads(e) for e in experiment_names)

    def get_reads(self, experiment_name):
        session = inspect(self).session
        r = session.query(Reads).get((experiment_name, self.gene_id, self.start, self.end))
        try:
            return r.reads
        except AttributeError:
            return 0

    def get_experiment(self):
        j = dict(self)
        j['intron_retention'] = self.get_intron_retention_type()
        return j

    def get_intron_retention_type(self):
        if self.intron_retention:
            for exon in filter(lambda e: not e.intron_retention, self.gene.exons):
                if self.start in exon:
                    return constants.IR_TYPE_START
                elif self.end in exon:
                    return constants.IR_TYPE_END

        return self.intron_retention


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

    @property
    def view_start(self):
        if self.start == -1:
            return self.end - 10
        return self.start

    @property
    def view_end(self):
        if self.end == -1:
            return self.start + 10
        return self.end

    def has_reads(self, experiment_name):
        return any(j.get_reads(experiment_name) > 0 for j in self.a5) or any(
            j.get_reads(experiment_name) > 0 for j in self.a3)

    def get_exon_type(self, experiment_name):
        if self.start == -1:
            return constants.EXON_TYPE_MISSING_START
        if self.end == -1:
            return constants.EXON_TYPE_MISSING_END

        has_reads = self.has_reads(experiment_name)

        if self.annotated and not has_reads:
            return constants.EXON_TYPE_DB

        if self.annotated and has_reads:
            return constants.EXON_TYPE_DB_RNASEQ

        if not self.annotated and has_reads:
            return constants.EXON_TYPE_RNASEQ

        return constants.EXON_TYPE_RNASEQ

    def get_experiment(self):
        exon = dict(self)
        exon['start'] = self.view_start
        exon['end'] = self.view_end
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

    def lsv_reference_exon(self, lsv_id):
        lsv_coords = lsv_id.split(':')[-1]
        if lsv_coords.startswith('-1'):
            coords_list = [-1, int(lsv_coords.split('-')[-1])]
        elif lsv_coords.endswith('-1'):
            coords_list = [int(lsv_coords.split('-')[0]), -1]
        else:
            coords_list = list(map(int, lsv_coords.split('-')))

        for exon in self.exons:
            if coords_list == [exon.start, exon.end]:
                return exon

    def lsv_exons(self, lsv, lsv_junctions=None):
        lsv_id = lsv.lsv_id
        if lsv_junctions is None:
            lsv_junctions = self.lsv_junctions(lsv)
        is_target = lsv_id.split(':')[-2] == 't'

        def find_exons():
            yield self.lsv_reference_exon(lsv_id)

            for junc in lsv_junctions:
                for exon in self.exons:
                    if is_target:
                        if self.strand == '+':
                            if junc.start in exon:
                                yield exon
                        else:
                            if junc.end in exon:
                                yield exon
                    else:
                        if self.strand == '+':
                            if junc.end in exon:
                                yield exon
                        else:
                            if junc.start in exon:
                                yield exon

        yield from sorted(find_exons(), key=lambda e: [e.start, e.end])

    def lsv_junctions(self, lsv):
        for start, end in lsv.junctions:
            for junc in self.junctions:
                if [start, end] == [junc.start, junc.end]:
                    yield junc

    def lsv_ucsc_coordinates(self, lsv):
        exons = tuple(self.lsv_exons(lsv))
        start_exon = sorted(e.start for e in exons if e.start != -1)[0]
        end_exon = sorted(e.end for e in exons if e.end != -1)[-1]
        return {'start': start_exon, 'end': end_exon}

    def get_experiment(self, experiment_names_list):
        gene = dict(self)
        gene['_id'] = self.id

        # todo: exons should NOT be sorted.
        exons = tuple(e.get_experiment() for e in self.exons)
        exons = sorted(exons, key=lambda e: (e['start'], e['end']))
        gene['exons'] = exons

        gene['junctions'] = tuple(j.get_experiment() for j in self.junctions)

        gene['start'] = self.start
        gene['end'] = self.end

        gene['exon_types'] = {}
        gene['reads'] = {}
        gene['junction_types'] = {}

        for experiment_names in experiment_names_list:
            combined_name = next((n for n in experiment_names if ' Combined' in n), '')
            experiment_names = experiment_names[experiment_names != combined_name]

            for exon in self.exons:
                for experiment_name in experiment_names:
                    try:
                        exon_start = gene['exon_types'][exon.view_start]
                    except KeyError:
                        gene['exon_types'][exon.view_start] = {}
                        exon_start = gene['exon_types'][exon.view_start]

                    try:
                        exon_end = exon_start[exon.view_end]
                    except KeyError:
                        exon_start[exon.view_end] = {}
                        exon_end = exon_start[exon.view_end]

                    exon_end[experiment_name] = exon.get_exon_type(experiment_name)

                if combined_name:
                    exon_end[combined_name] = min(exon_end[x] for x in experiment_names)

            for junc in self.junctions:

                for experiment_name in experiment_names:
                    try:
                        reads_start = gene['reads'][junc.start]
                    except KeyError:
                        gene['reads'][junc.start] = {}
                        reads_start = gene['reads'][junc.start]

                    try:
                        reads_end = reads_start[junc.end]
                    except KeyError:
                        reads_start[junc.end] = {}
                        reads_end = reads_start[junc.end]

                    reads_end[experiment_name] = junc.get_reads(experiment_name)

                try:
                    types_start = gene['junction_types'][junc.start]
                except KeyError:
                    gene['junction_types'][junc.start] = {}
                    types_start = gene['junction_types'][junc.start]

                try:
                    types_end = types_start[junc.end]
                except KeyError:
                    types_start[junc.end] = {}
                    types_end = types_start[junc.end]

                for name, junc_type in junc.get_junction_types(experiment_names):
                    types_end[name] = junc_type

                if combined_name:
                    reads = gene['reads'][junc.start][junc.end]
                    junction_types = gene['junction_types'][junc.start][junc.end]
                    reads[combined_name] = sum(reads[x] for x in experiment_names)
                    junction_types[combined_name] = min(junction_types[x] for x in experiment_names)

        return gene
