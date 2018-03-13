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

    def get_experiment(self):
        j = dict(self)
        j['intron_retention'] = self.get_intron_retention_type()
        j['color'] = self.color()
        del j['gene_id']
        return j

    def get_intron_retention_type(self):
        if self.intron_retention:
            for exon in filter(lambda e: not e.intron_retention, self.gene.exons):
                if self.start in exon:
                    return constants.IR_TYPE_START
                elif self.end in exon:
                    return constants.IR_TYPE_END

        return self.intron_retention


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

    def get_experiment(self):
        exon = dict(self)
        exon['start'] = self.view_start
        exon['end'] = self.view_end
        if self.start == -1:
            exon['half_exon'] = 'start'
        elif self.end == -1:
            exon['half_exon'] = 'end'
        exon['color'] = self.color()
        del exon['gene_id']
        return exon


class ViewGene(Gene):
    junctions = relationship('ViewJunction')
    exons = relationship('ViewExon')

    @staticmethod
    def convert_lsv_to_coords(lsv_id):
        for coord in lsv_id.split(':')[-1].split('-'):
            if coord == 'nan':
                yield -1
            else:
                yield int(coord)

    def lsv_reference_exon(self, lsv_id):
        coords_list = list(self.convert_lsv_to_coords(lsv_id))
        return next(exon for exon in self.exons if coords_list == [exon.start, exon.end])

    def lsv_exons(self, lsv):
        exons = {exon for start, end in lsv.junctions for exon in self.exons if start in exon or end in exon}
        yield from sorted(exons, key=lambda e: [e.start, e.end])

    def lsv_junctions(self, lsv):
        for start, end in lsv.junctions:
            yield next(junc for junc in self.junctions if [start, end] == [junc.start, junc.end])

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
        gene['reads'] = {}

        for experiment_names in experiment_names_list:
            combined_name = next((n for n in experiment_names if ' Combined' in n), '')
            experiment_names = experiment_names[experiment_names != combined_name]

            for name in experiment_names:
                gene['reads'][name] = {}
                if combined_name:
                    gene['reads'][combined_name] = {}

            for junc in self.junctions:
                view_junc = junc
                for experiment_name in experiment_names:
                    try:
                        gene['reads'][experiment_name][junc.start][junc.end] = view_junc.reads(experiment_name)
                    except KeyError:
                        gene['reads'][experiment_name][junc.start] = {junc.end: view_junc.reads(experiment_name)}

                if combined_name:
                    summed_reads = sum(gene['reads'][n][junc.start][junc.end] for n in experiment_names)
                    try:
                        gene['reads'][combined_name][junc.start][junc.end] = summed_reads
                    except KeyError:
                        gene['reads'][combined_name][junc.start] = {junc.end: summed_reads}

        return gene
