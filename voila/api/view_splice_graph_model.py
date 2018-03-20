from sqlalchemy import inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from voila import constants
from voila.api.splice_graph_model import Gene, Exon, Junction, Reads

Base = declarative_base()


class ViewJunction(Junction):
    @property
    def id(self):
        return '{}:{}-{}'.format(self.gene.id, self.start, self.end)

    def color(self):
        if self.annotated:
            if self.has_reads:
                return 'red'
            else:
                return 'grey'
        else:
            return 'green'

    def reads(self, experiment_name):
        session = inspect(self).session
        r = session.query(Reads).get((experiment_name, self.gene_id, self.start, self.end))
        try:
            return r.reads
        except AttributeError:
            return 0

    def view_intron_retention_type(self):
        if self.intron_retention:
            for exon in filter(lambda e: not e.intron_retention, self.gene.exons):
                if self.start in exon:
                    return constants.IR_TYPE_START
                elif self.end in exon:
                    return constants.IR_TYPE_END

        return self.intron_retention

    def __iter__(self):
        yield 'start', self.start
        yield 'end', self.end
        yield 'has_reads', self.has_reads
        yield 'intron_retention', self.view_intron_retention_type()
        yield 'annotated', self.annotated
        yield 'color', self.color()


class ViewExon(Exon):
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

    def color(self):
        if self.annotated:
            if any(j.has_reads for j in self.a5) or any(j.has_reads for j in self.a3):
                return 'grey'
            else:
                return ''

        else:
            return 'green'

    def __iter__(self):
        yield 'start', self.view_start
        yield 'end', self.view_end
        if self.start == -1:
            yield 'half_exon', 'start'
        elif self.end == -1:
            yield 'half_exon', 'end'
        yield 'intron_retention', self.intron_retention
        yield 'annotated', self.annotated
        yield 'color', self.color()


class ViewGene(Gene):
    junctions = relationship('ViewJunction')
    exons = relationship('ViewExon')

    def __iter__(self):
        yield 'name', self.name
        yield 'strand', self.strand
        yield 'chromosome', self.chromosome
        yield '_id', self.id
        yield 'id', self.id
        yield 'start', self.start
        yield 'end', self.end

    @staticmethod
    def convert_lsv_to_coords(lsv_id):
        for coord in lsv_id.split(':')[-1].split('-'):
            if coord == 'nan':
                yield -1
            else:
                yield int(coord)

    def lsv_exons(self, lsv):
        exons = {exon for start, end in lsv.junctions for exon in self.exons if start in exon or end in exon}
        yield from sorted(exons, key=lambda e: [e.start, e.end])

    def lsv_junctions(self, lsv):
        for start, end in lsv.junctions:
            yield next(junc for junc in self.junctions if [start, end] == [junc.start, junc.end])

    def lsv_ucsc_coordinates(self, lsv):
        exons = list(self.lsv_exons(lsv))
        start_exon = sorted(e.start for e in exons if e.start != -1)[0]
        end_exon = sorted(e.end for e in exons if e.end != -1)[-1]
        return {'start': start_exon, 'end': end_exon}

    def get_experiment(self, experiment_names_list):
        reads = {}

        for experiment_names in experiment_names_list:
            combined_name = next((n for n in experiment_names if ' Combined' in n), '')
            experiment_names = experiment_names[experiment_names != combined_name]

            for name in experiment_names:
                reads[name] = {}
                if combined_name:
                    reads[combined_name] = {}

            for junc in self.junctions:
                for experiment_name in experiment_names:
                    try:
                        reads[experiment_name][junc.start][junc.end] = junc.reads(experiment_name)
                    except KeyError:
                        reads[experiment_name][junc.start] = {junc.end: junc.reads(experiment_name)}

                if combined_name:
                    summed_reads = sum(reads[n][junc.start][junc.end] for n in experiment_names)
                    try:
                        reads[combined_name][junc.start][junc.end] = summed_reads
                    except KeyError:
                        reads[combined_name][junc.start] = {junc.end: summed_reads}

        gene = dict(self)
        gene['exons'] = [dict(e) for e in self.exons]
        gene['junctions'] = [dict(j) for j in self.junctions]
        gene['reads'] = reads

        return gene
