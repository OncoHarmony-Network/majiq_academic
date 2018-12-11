from operator import itemgetter

from voila.api import SpliceGraph
from voila.config import ViewConfig


class ViewSpliceGraph(SpliceGraph):
    def __init__(self):
        config = ViewConfig()
        splice_graph_file = config.splice_graph_file
        super().__init__(splice_graph_file)

    @property
    def metadata(self):
        query = self.conn.execute('SELECT name FROM experiment')
        experiment_names = [x for x, in query.fetchall()]
        return {'experiment_names': [experiment_names], 'group_names': ['Splice Graph']}

    @staticmethod
    def exon_start(exon):
        if exon['start'] == -1:
            return exon['end'] - 10
        return exon['start']

    @staticmethod
    def exon_end(exon):
        if exon['end'] == -1:
            return exon['start'] + 10
        return exon['end']

    @staticmethod
    def exon_annot_start(exon):
        if exon['annotated_start'] == -1:
            return exon['annotated_end'] - 10
        return exon['annotated_start']

    @staticmethod
    def exon_annot_end(exon):
        if exon['annotated_end'] == -1:
            return exon['annotated_start'] + 10
        return exon['annotated_end']

    @property
    def gene_ids(self):
        query = self.conn.execute('SELECT id FROM gene')
        return [x for x, in query.fetchall()]

    def view_gene(self, gene_id):
        gene = self.gene(gene_id)
        yield from gene.items()
        yield 'id', gene_id
        yield 'start', self.gene_start(gene_id)
        yield 'end', self.gene_end(gene_id)

    def gene_start(self, gene_id):
        return sorted(self.exon_start(e) for e in self.exons(gene_id))[0]

    def gene_end(self, gene):
        return sorted((self.exon_end(e) for e in self.exons(gene)), reverse=True)[0]

    def view_exon(self, exon):
        yield 'start', self.exon_start(exon)
        yield 'end', self.exon_end(exon)
        if exon['start'] == -1:
            yield 'half_exon', 'start'
        elif exon['end'] == -1:
            yield 'half_exon', 'end'
        yield 'annotated', exon['annotated']
        yield 'color', self.exon_color(exon)
        yield 'annotated_start', self.exon_annot_start(exon)
        yield 'annotated_end', self.exon_annot_end(exon)

    def view_exons(self, gene_id):
        for exon in self.exons(gene_id):
            yield self.view_exon(exon)

    def exon_has_reads(self, exon):
        query = self.conn.execute('''
                        SELECT has_reads FROM junction
                        WHERE 
                        (gene_id=? AND has_reads=1)
                        AND 
                        (
                          ({0}=-1 AND start={1})
                          OR 
                          ({1}=-1 AND end={0})
                          OR
                          (-1 NOT IN ({0},{1}) AND start BETWEEN {0} AND {1})
                          OR 
                          (-1 NOT IN ({0},{1}) AND end BETWEEN {0} AND {1})
                        )
                        LIMIT 1
                      '''.format(exon['start'], exon['end']), (exon['gene_id'],))
        return query.fetchone()

    def exon_color(self, exon):
        if exon['annotated']:
            if self.exon_has_reads(exon):
                return 'grey'
            else:
                return ''
        else:
            return 'green'

    def view_junctions(self, gene):
        for junc in self.junctions(gene):
            yield self.view_junction(junc)

    def view_junction(self, junction):
        yield from junction.items()
        yield 'color', self.junction_color(junction)

    def junction_color(self, junction):
        if junction['annotated']:
            if junction['has_reads']:
                return 'red'
            else:
                return 'grey'
        else:
            return 'green'

    def view_intron_retentions(self, gene):
        for ir in self.intron_retentions(gene):
            yield self.view_intron_retention(ir)

    def view_intron_retention(self, ir):
        yield from ir.items()
        yield 'color', self.ir_color(ir)

    def ir_color(self, ir):
        if ir['annotated']:
            if ir['has_reads']:
                return 'red'
            else:
                return 'grey'
        else:
            return 'green'

    def annotated_junctions(self, gene_id, lsv_junctions):
        for junc in lsv_junctions:
            junc = tuple(map(int, junc))
            query = self.conn.execute('''
                                SELECT annotated FROM junction
                                WHERE gene_id=?
                                AND start=? 
                                AND end=?
                                ''', (gene_id, junc[0], junc[1]))

            fetch = query.fetchone()
            if fetch:
                yield fetch[0]

    def lsv_exons(self, gene_id, lsv_junctions):
        rtn_set = set()
        for junc in lsv_junctions:
            junc = tuple(map(int, junc))
            query = self.conn.execute('''
                                        SELECT start, end FROM exon
                                        WHERE gene_id=? 
                                        AND
                                        (
                                        (start=-1 AND end=?)
                                        OR 
                                        (end=-1 AND start=?)
                                        OR
                                        (start!=-1 AND end!=-1 AND ? BETWEEN start AND end)
                                        OR 
                                        (start!=-1 AND end!=-1 AND ? BETWEEN start and end)
                                        )  
                                        ''', (gene_id, junc[0], junc[1], junc[0], junc[1]))

            for x in query.fetchall():
                rtn_set.add(x)

        return list(sorted(rtn_set))

    def lsv_introns(self, gene, lsv_exons):
        exons_ends = list(e for s, e in lsv_exons)
        for ir in self.intron_retentions(gene):
            if ir.start in exons_ends:
                yield ir

    def gene_experiment(self, gene_id, experiment_names_list):
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

            for junc in self.junctions(gene_id):
                junc_start, junc_end = itemgetter('start', 'end')(junc)

                for r in self.junction_reads_exp(junc, experiment_names):
                    reads = r['reads']
                    exp_name = r['experiment_name']
                    try:
                        junc_reads[exp_name][junc_start][junc_end] = reads
                    except KeyError:
                        junc_reads[exp_name][junc_start] = {junc_end: reads}

                    if combined_name:
                        def get_junc_reads():
                            for n in experiment_names:
                                try:
                                    yield junc_reads[n][junc_start][junc_end]
                                except KeyError:
                                    pass

                        summed_reads = sum(get_junc_reads())

                        try:
                            junc_reads[combined_name][junc_start][junc_end] = summed_reads
                        except KeyError:
                            junc_reads[combined_name][junc_start] = {junc_end: summed_reads}

            for ir in self.intron_retentions(gene_id):

                ir_start, ir_end = itemgetter('start', 'end')(ir)

                for r in self.intron_retention_reads_exp(ir, experiment_names):

                    reads, exp_name = itemgetter('reads', 'experiment_name')(r)

                    try:
                        ir_reads[exp_name][ir_start][ir_end] = reads
                    except KeyError:
                        ir_reads[exp_name][ir_start] = {ir_end: reads}

                if combined_name:
                    def get_ir_reads():
                        for n in experiment_names:
                            try:
                                yield ir_reads[n][ir_start][ir_end]
                            except KeyError:
                                pass

                    summed_reads = sum(get_ir_reads())

                    try:
                        ir_reads[combined_name][ir_start][ir_end] = summed_reads
                    except KeyError:
                        ir_reads[combined_name][ir_start] = {ir_end: summed_reads}

        gene_dict = dict(self.view_gene(gene_id))
        gene_dict['exons'] = tuple(dict(e) for e in self.view_exons(gene_id))
        gene_dict['junctions'] = tuple(dict(j) for j in self.view_junctions(gene_id))
        gene_dict['intron_retention'] = tuple(dict(ir) for ir in self.view_intron_retentions(gene_id))
        gene_dict['junction_reads'] = junc_reads
        gene_dict['intron_retention_reads'] = ir_reads
        gene_dict['genome'] = self.genome
        gene_dict['alt_starts'] = tuple(list(a.values())[0] for a in self.alt_starts(gene_id))
        gene_dict['alt_ends'] = tuple(list(a.values())[0] for a in self.alt_ends(gene_id))

        return gene_dict
