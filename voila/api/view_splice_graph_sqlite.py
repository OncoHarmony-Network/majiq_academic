from voila.api import SpliceGraph


class ViewSpliceGraph(SpliceGraph):
    def __init__(self, args):
        super().__init__(args.splice_graph)

    def gene(self, gene_id):
        print(super().gene(gene_id))

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

        gene_dict = gene._asdict()
        gene_dict['exons'] = [e._asdict() for e in self.exons(gene)]
        gene_dict['junctions'] = [j._asdict() for j in self.junctions(gene)]
        gene_dict['intron_retention'] = [ir._asdict() for ir in self.intron_retentions(gene)]
        gene_dict['junction_reads'] = junc_reads
        gene_dict['intron_retention_reads'] = ir_reads
        gene_dict['genome'] = self.genome

        return gene_dict
