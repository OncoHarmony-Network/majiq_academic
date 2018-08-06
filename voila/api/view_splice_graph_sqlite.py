from voila.api import SpliceGraph


class ViewSpliceGraph(SpliceGraph):
    def __init__(self, args):
        super().__init__(args.splice_graph)

    @property
    def metadata(self):
        query = self.conn.execute('''
                            SELECT name FROM experiment
                            ''')
        experiment_names = [x for x, in query.fetchall()]
        return {'experiment_names': [experiment_names], 'group_names': ['Splice Graph'], '_id': 'metadata'}

    def view_gene(self, gene):
        yield from gene._asdict().items()
        yield '_id', gene.id
        yield 'start', self.gene_start(gene)
        yield 'end', self.gene_end(gene)

    def gene_start(self, gene):
        return sorted(self.exon_start(e) for e in self.exons(gene))[0]

    def gene_end(self, gene):
        return sorted((self.exon_end(e) for e in self.exons(gene)), reverse=True)[0]

    @staticmethod
    def exon_start(exon):
        if exon.start == -1:
            return exon.end - 10
        return exon.start

    @staticmethod
    def exon_end(exon):
        if exon.end == -1:
            return exon.start + 10
        return exon.end

    def view_exon(self, exon):
        yield 'start', self.exon_start(exon)
        yield 'end', self.exon_end(exon)
        if exon.start == -1:
            yield 'half_exon', 'start'
        elif exon.end == -1:
            yield 'half_exon', 'end'
        yield 'annotated', exon.annotated
        yield 'color', self.exon_color(exon)

    def view_exons(self, gene):
        for exon in self.exons(gene):
            yield dict(self.view_exon(exon))

    def exon_has_reads(self, exon):
        query = self.conn.execute('''
                        SELECT has_reads FROM junction
                        WHERE 
                        (has_reads=1)
                        AND 
                        (
                          ({0}=-1 AND start={1})
                          OR 
                          ({1}=-1 AND end={0})
                          OR
                          (-1 NOT IN ({0},{1}) AND start BETWEEN {0} AND {1})
                          OR 
                          (-1 NOT IN ({0},{1}) AND END BETWEEN {0} AND {1})
                        )
                        LIMIT 1
                      '''.format(exon.start, exon.end))
        return query.fetchone()

    def exon_color(self, exon):
        if exon.annotated:
            if self.exon_has_reads(exon):
                return 'grey'
            else:
                return ''
        else:
            return 'green'

    def view_junction(self, junction):
        yield from junction._asdict().items()
        yield 'color', self.junction_color(junction)

    def view_junctions(self, gene):
        for junc in self.junctions(gene):
            yield dict(self.view_junction(junc))

    def junction_color(self, junction):
        if junction.annotated:
            if junction.has_reads:
                return 'red'
            else:
                return 'grey'
        else:
            return 'green'

    def view_intron_retention(self, ir):
        yield from ir._asdict().items()
        yield 'color', self.ir_color(ir)

    def view_intron_retentions(self, gene):
        for ir in self.intron_retentions(gene):
            yield dict(self.view_intron_retention(ir))

    def ir_color(self, ir):
        if ir.annotated:
            if ir.has_reads:
                return 'red'
            else:
                return 'grey'
        else:
            return 'green'

    def gene_experiment(self, gene, experiment_names_list):
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

            for junc in self.junctions(gene):
                reads = (r for r in self.junction_reads(junc) if r.experiment_name in experiment_names)
                for r in reads:
                    try:
                        junc_reads[r.experiment_name][junc.start][junc.end] = r.reads
                    except KeyError:
                        junc_reads[r.experiment_name][junc.start] = {junc.end: r.reads}

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

            for ir in self.intron_retentions(gene):
                reads = (r for r in self.intron_retention_reads(ir) if r.experiment_name in experiment_names)
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

        gene_dict = dict(self.view_gene(gene))
        gene_dict['exons'] = tuple(self.view_exons(gene))
        gene_dict['junctions'] = tuple(self.view_junctions(gene))
        gene_dict['intron_retention'] = tuple(self.view_intron_retentions(gene))
        gene_dict['junction_reads'] = junc_reads
        gene_dict['intron_retention_reads'] = ir_reads
        gene_dict['genome'] = self.genome

        return gene_dict
