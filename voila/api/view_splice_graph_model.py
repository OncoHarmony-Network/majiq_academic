from sqlalchemy import inspect

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from voila.api.splice_graph_model import Gene, Exon, Junction, JunctionReads, IntronRetention, Genome

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

    def __iter__(self):
        yield 'start', self.start
        yield 'end', self.end
        yield 'has_reads', self.has_reads
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
        yield 'annotated', self.annotated
        yield 'color', self.color()


class ViewIntronRetention(IntronRetention):
    def __iter__(self):
        yield 'start', self.start
        yield 'end', self.end
        yield 'annotated', self.annotated
        yield 'color', self.color()

    def color(self):
        if self.annotated:
            if self.has_reads:
                return 'red'
            else:
                return 'grey'
        else:
            return 'green'


class ViewGene(Gene):
    junctions = relationship('ViewJunction')
    exons = relationship('ViewExon')
    intron_retentions = relationship('ViewIntronRetention')

    def __iter__(self):
        yield 'name', self.name
        yield 'strand', self.strand
        yield 'chromosome', self.chromosome
        yield 'id', self.id
        yield '_id', self.id
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
        junc_reads = {}
        ir_reads = {}

        for experiment_names in experiment_names_list:
            combined_name = next((n for n in experiment_names if ' Combined' in n), '')
            experiment_names = [e for e in experiment_names if e != combined_name]

            for name in experiment_names:
                junc_reads[name] = {}
                ir_reads[name] = {}
                if combined_name:
                    junc_reads[combined_name] = {}
                    ir_reads[combined_name] = {}

            for junc in self.junctions:
                reads = (r for r in junc.reads if r.experiment_name in experiment_names)
                for r in reads:
                    try:
                        junc_reads[r.experiment_name][r.junction_start][r.junction_end] = r.reads
                    except KeyError:
                        junc_reads[r.experiment_name][r.junction_start] = {r.junction_end: r.reads}

                    if combined_name:
                        def get_junc_reads():
                            for n in experiment_names:
                                try:
                                    yield junc_reads[n][junc.start][junc.end]
                                except KeyError:
                                    pass

                        summed_reads = sum(get_junc_reads())

                        try:
                            junc_reads[combined_name][junc.start][junc.end] = summed_reads
                        except KeyError:
                            junc_reads[combined_name][junc.start] = {junc.end: summed_reads}

            for ir in self.intron_retentions:
                reads = (r for r in ir.reads if r.experiment_name in experiment_names)
                for r in reads:
                    try:
                        ir_reads[r.experiment_name][ir.start][ir.end] = r.reads
                    except KeyError:
                        ir_reads[r.experiment_name][ir.start] = {ir.end: r.reads}

                if combined_name:
                    def get_ir_reads():
                        for n in experiment_names:
                            try:
                                yield ir_reads[n][ir.start][ir.end]
                            except KeyError:
                                pass

                    summed_reads = sum(get_ir_reads())

                    try:
                        ir_reads[combined_name][ir.start][ir.end] = summed_reads
                    except KeyError:
                        ir_reads[combined_name][ir.start] = {ir.end: summed_reads}

        gene = dict(self)
        gene['exons'] = [dict(e) for e in self.exons]
        gene['junctions'] = [dict(j) for j in self.junctions]
        gene['intron_retention'] = [dict(ir) for ir in self.intron_retentions]
        gene['junction_reads'] = junc_reads
        gene['intron_retention_reads'] = ir_reads
        gene['genome'] = inspect(self).session.query(Genome).one().name

        return gene
