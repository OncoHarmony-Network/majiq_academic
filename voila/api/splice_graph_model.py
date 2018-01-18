from sqlalchemy import Integer, String, Column, Boolean, ForeignKey, ForeignKeyConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from voila import constants

Base = declarative_base()


def coord_in_exon(coord, exon):
    if -1 in (exon.start, exon.end):
        return coord in (exon.start, exon.end)
    else:
        return exon.start <= coord <= exon.end


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
        junction['intron_retention'] = self.get_intron_retention_type()

        return junction

    def get_intron_retention_type(self):
        if self.intron_retention:
            for exon in (e for e in self.gene.exons if not e.intron_retention):
                if self.start == exon.end:
                    return constants.IR_TYPE_START
                elif self.end == exon.start:
                    return constants.IR_TYPE_END

        return self.intron_retention


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
    def a3(self):
        def a3_filter(js):
            if self.end == -1:
                return

            for j in js:
                if coord_in_exon(j.start, self):
                    yield j

        return list(a3_filter(self.gene.junctions))

    @property
    def a5(self):
        def a5_filter(js):
            if self.start == -1:
                return

            for j in js:
                if coord_in_exon(j.end, self):
                    yield j

        return list(a5_filter(self.gene.junctions))

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

    def lsv_exons(self, lsv_id, lsv_junctions=None):
        is_target = lsv_id.split(':')[-2] == 't'
        if lsv_junctions is None:
            lsv_junctions = self.lsv_junctions(lsv_id)

        ref_exon = self.lsv_reference_exon(lsv_id)

        exons = {ref_exon}
        for junc in lsv_junctions:
            for exon in self.exons:
                if is_target:
                    if self.strand == '+':
                        if coord_in_exon(junc.start, exon):
                            exons.add(exon)
                    else:
                        if coord_in_exon(junc.end, exon):
                            exons.add(exon)
                else:
                    if self.strand == '+':
                        if coord_in_exon(junc.end, exon):
                            exons.add(exon)
                    else:
                        if coord_in_exon(junc.start, exon):
                            exons.add(exon)

        assert len(set('{}-{}'.format(exon.start, exon.end) for exon in exons)) == len(exons)

        return exons

    def lsv_junctions(self, lsv_id):
        is_target = lsv_id.split(':')[-2] == 't'
        exon = self.lsv_reference_exon(lsv_id)
        if is_target:
            if self.strand == '+':
                return exon.a5
            else:
                return exon.a3
        else:
            if self.strand == '+':
                return exon.a3
            else:
                return exon.a5

    def lsv_coordinates(self, lsv_id):
        d = lsv_id.split(self.id)[1].split(':')[1:]
        try:
            coordinates = list(map(int, d[1].split('-')))
        except ValueError:
            coordinates = [-1, int(d[1].split('-')[-1])]

        start = coordinates[0]
        end = coordinates[1]
        is_target = d[0] == 't'

        # print('coords', coordinates)

        if self.strand == '-':
            return {'start': 0, 'end': 0}

        for exon in self.exons:
            if exon.start == start and exon.end == end:
                # print('start', exon.start)
                # print('end', exon.end)
                # print('is_target', is_target)

                if is_target:
                    # print('using a5')
                    juncs = exon.a5
                    j = sorted(juncs, key=lambda junc: junc.start)[0]
                    for start_exon in self.exons:
                        # if start_exon.start <= j.start <= start_exon.end:
                        if coord_in_exon(j.start, start_exon):
                            return {'start': start_exon.view_start, 'end': exon.view_end}
                else:
                    # print('using a3')
                    juncs = exon.a3
                    j = sorted(juncs, key=lambda junc: junc.end)[-1]

                    for end_exon in self.exons:
                        # if end_exon.start <= j.end <= end_exon.end:
                        if coord_in_exon(j.end, end_exon):
                            return {'start': exon.view_start, 'end': end_exon.view_end}

        return {'start': 0, 'end': 0}

    def get_experiment(self, experiment_names):
        gene = dict(self)

        gene['_id'] = self.id

        # todo: exons should NOT be sorted.
        exons = tuple(e.get_experiment() for e in self.exons)
        exons = sorted(exons, key=lambda e: (e['start'], e['end']))
        gene['exons'] = exons

        gene['junctions'] = tuple(dict(j) for j in self.junctions)

        gene['start'] = self.start
        gene['end'] = self.end

        gene['exon_types'] = {exon.view_start: {} for exon in self.exons}
        for exon in self.exons:
            gene['exon_types'][exon.view_start][exon.view_end] = {
                experiment_name: exon.get_exon_type(experiment_name) for experiment_name in experiment_names
            }

        gene['reads'] = {junc.start: {} for junc in self.junctions}
        gene['junction_types'] = {junc.start: {} for junc in self.junctions}
        for junc in self.junctions:
            gene['reads'][junc.start][junc.end] = {
                experiment_name: junc.get_reads(experiment_name) for experiment_name in experiment_names
            }
            gene['junction_types'][junc.start][junc.end] = {
                experiment_name: junc.get_junction_type(experiment_name) for experiment_name in experiment_names
            }

        return gene
