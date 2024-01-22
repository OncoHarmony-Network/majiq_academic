import sqlite3
from operator import itemgetter

from rna_voila.api import SpliceGraph
from rna_voila.api.splice_graph import transcript_exon_fieldnames
from rna_voila.config import ViewConfig
from rna_voila.api.splice_graph_lr import combined_colors
from statistics import median, StatisticsError
from math import ceil


class ViewSpliceGraph(SpliceGraph):
    def __init__(self, omit_simplified=False, splice_graph_file=None):
        """
        Wrapper class to splice graph api used to generate voila output.
        """
        self.omit_simplified = omit_simplified
        if not splice_graph_file:
            config = ViewConfig()
            splice_graph_file = config.splice_graph_file
        super().__init__(splice_graph_file)

    @staticmethod
    def exon_start(exon):
        """
        Get start of exon. If no start, then use end to calculate start.
        :param exon: exon dictionary from splice graph
        :return: integer
        """

        if exon['start'] == -1:
            return exon['end'] - 10
        return exon['start']

    @staticmethod
    def exon_end(exon):
        """
        Get end of exon. If no end, then use start to calculate end.
        :param exon: exon dictionary from splice graph
        :return: integer
        """

        if exon['end'] == -1:
            return exon['start'] + 10
        return exon['end']

    @staticmethod
    def exon_annot_start(exon):
        """
        Get annotated start. When exon is denovo/missing start, then this value might be -1.
        :param exon: exon dictionary from splice graph
        :return: integer
        """

        if exon['annotated_start'] == -1:
            return exon['annotated_end'] - 10
        return exon['annotated_start']

    @staticmethod
    def exon_annot_end(exon):
        """
        Get annotated end. When exon is denovo/missing end, then this value might be -1.
        :param exon: exon dictionary from splice graph
        :return: integer
        """

        if exon['annotated_end'] == -1:
            return exon['annotated_start'] + 10
        return exon['annotated_end']

    @property
    def gene_ids(self):
        """
        List of all gene ids in splice graph.
        :return: list
        """

        query = self.conn.execute('SELECT id FROM gene')
        return [x for x, in query.fetchall()]

    @property
    def gene_ids2gene_names(self):
        """
        Dict of all gene_ids to gene_names
        :return: list
        """

        query = self.conn.execute('SELECT id, name FROM gene')
        return {x[0]: x[1] for x in query.fetchall()}

    def view_gene(self, gene_id):
        """
        Add information to gene dictionary that is used by javascript splice graph.
        :param gene_id: gene id
        :return: generator key/value
        """

        gene = self.gene(gene_id)
        yield from gene.items()
        yield 'id', gene_id
        yield 'start', self.gene_start(gene_id)
        yield 'end', self.gene_end(gene_id)

    def gene_start(self, gene_id):
        """
        Find gene start from exon coords.
        :param gene_id: gene id
        :return: integer
        """

        return sorted(self.exon_start(e) for e in self.exons(gene_id))[0]

    def gene_end(self, gene):
        """
        Find gene end from exon coords.
        :param gene_id: gene id
        :return: integer
        """

        return sorted((self.exon_end(e) for e in self.exons(gene)), reverse=True)[0]

    def view_exon(self, exon):
        """
        Add information to exon dictionary that is used by javascript splice graph.
        :param exon: exon dictionary from splice graph file.
        :return: generator key/value
        """

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
        """
        Get all view exons.
        :param gene_id: gene id
        :return: generator of view exons
        """

        for exon in self.exons(gene_id):
            yield self.view_exon(exon)

    def exon_has_reads(self, exon):
        """
        Does this exon have at least one junction that has reads.
        :param exon: exon dictionary from splice graph
        :return: boolean
        """

        found = self.conn.execute('''
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
                        {2}
                        LIMIT 1
                      '''.format(exon['start'], exon['end'],
                                (" AND is_simplified = 0" if self.omit_simplified else '')),
                                (exon['gene_id'],)).fetchone()
        if not found:
            # if no reads in the junction table, also check the intron retention table
            found = self.conn.execute('''
                SELECT has_reads FROM intron_retention
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
                {2}
                LIMIT 1
              '''.format(exon['start'] - 1, exon['end'] + 1,
                         (" AND is_simplified = 0" if self.omit_simplified else '')),
                  (exon['gene_id'],)).fetchone()

        return found

    def exon_color(self, exon):
        """
        Get color of exon for javascript splice graph.
        :param exon: exon dictionary from splice graph
        :return: string
        """

        if exon['annotated']:
            if self.exon_has_reads(exon):
                return combined_colors['ao']
            else:
                return ''
        else:
            return combined_colors['s']

    def view_junctions(self, gene):
        """
        Get all view junctions.
        :param gene: gene dictionary
        :return: generator of all junctions
        """

        for junc in self.junctions(gene, omit_simplified=self.omit_simplified):
            yield self.view_junction(junc)

    def view_junction(self, junction):
        """
        Add information to juction dictionary for javascript splice graph.
        :param junction: junction dictionary from splice graph file.
        :return: generator key/value
        """

        yield from junction.items()
        yield 'color', self.junction_color(junction)

    def junction_color(self, junction):
        """
        Find color for junction to use with javascript splice graph.
        :param junction: junction dictionary from splice graph file.
        :return: string
        """

        if junction['annotated']:
            if junction['has_reads']:
                return combined_colors['sa']
            else:
                return combined_colors['ao']
        else:
            return combined_colors['s']

    def view_intron_retentions(self, gene):
        """
        Get all view irs.
        :param gene: gene dictionary from splice graph.
        :return: generator
        """

        for ir in self.intron_retentions(gene, omit_simplified=self.omit_simplified):
            yield self.view_intron_retention(ir)

    def view_intron_retention(self, ir):
        """
        Add information to the ir dictionary from splice graph file.
        :param ir: ir dictionary from splice graph file.
        :return: generator key/value
        """

        yield from ir.items()
        yield 'color', self.ir_color(ir)

    def ir_color(self, ir):
        """
        Find color for intron retention.
        :param ir: ir dictionary from splice graph file.
        :return: string
        """

        if ir['annotated']:
            if ir['has_reads']:
                return combined_colors['sa']
            else:
                return combined_colors['ao']
        else:
            return combined_colors['s']

    def annotated_junctions(self, gene_id, lsv_junctions):
        """
        List of junction which are annotated in the db.
        :param gene_id: gene id
        :param lsv_junctions: list of juctions for an LSV.
        :return: generator
        """

        for junc in lsv_junctions:
            junc = tuple(map(int, junc))
            junc_query = self.conn.execute('''
                                SELECT annotated FROM junction
                                WHERE gene_id=?
                                AND start=? 
                                AND end=?
                                ''', (gene_id, junc[0], junc[1]))
            junc_res = junc_query.fetchone()
            if junc_res:
                yield junc_res[0]
            else:
                intron_query = self.conn.execute('''
                                                SELECT annotated FROM intron_retention
                                                WHERE gene_id=?
                                                AND start=? 
                                                AND end=?
                                                ''', (gene_id, junc[0], junc[1]))
                intron_res = intron_query.fetchone()
                if intron_res:
                    yield intron_res[0]

    def lsv_reads(self, gene_id, lsv_junctions, experiment_names, has_ir, dbg=False):
        """
        List of junction which are annotated in the db.
        :param gene_id: gene id
        :param lsv_junctions: list of juctions for an LSV.
        :return: generator
        """

        exps = {exp: [[], []] for exp in experiment_names}
        if has_ir:
            ir_junction = lsv_junctions[-1]
            lsv_junctions = lsv_junctions[:-1]
        else:
            ir_junction = None



        for junc in lsv_junctions:
            junc = tuple(map(int, junc))
            found_in_exps = set(exps.keys())
            junc_query = self.conn.execute(f'''
                                SELECT reads, experiment_name FROM junction_reads
                                WHERE junction_gene_id=?
                                AND junction_start=? 
                                AND junction_end=?
                                AND experiment_name in ({','.join(['?']*len(found_in_exps))})
                                ''', (gene_id, str(junc[0]), str(junc[1]), *found_in_exps))
            junc_res = junc_query.fetchall()
            for junc in junc_res:
                reads, exp = junc
                found_in_exps.remove(exp)
                exps[exp][0].append(reads)
            for exp in found_in_exps:
                exps[exp][0].append(0)

        if has_ir:
            found_in_exps = set(exps.keys())
            intron_query = self.conn.execute(f'''
                                            SELECT reads, experiment_name FROM intron_retention_reads
                                            WHERE intron_retention_gene_id=?
                                            AND intron_retention_start=? 
                                            AND intron_retention_end=?
                                            AND experiment_name in ({','.join(['?']*len(found_in_exps))})
                                            ''', (gene_id, str(ir_junction[0]), str(ir_junction[1]), *found_in_exps))

            intron_res = intron_query.fetchall()
            for intron in intron_res:
                reads, exp = intron

                found_in_exps.remove(exp)
                exps[exp][1].append(reads)
            for exp in found_in_exps:
                exps[exp][1].append(0)

        return exps


    def lsv_exons(self, gene_id, lsv_junctions):
        """
        Get exons for an LSV.
        :param gene_id: gene id
        :param lsv_junctions: list of lsv junctions
        :return: list
        """

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
        """
        Get ir for an LSV.
        :param gene: gene dictionary from splice graph file.
        :param lsv_exons: list of lsv exons
        :return: generator
        """
        exons_ends = list(e for s, e in lsv_exons)
        for ir in self.intron_retentions(gene, omit_simplified=self.omit_simplified):
            if ir.start in exons_ends:
                yield ir

    def _add_presence_key(self, items):
        for item in items:
            if not item['has_reads']:
                item['presence'] = 'ao'
            elif item['annotated']:
                item['presence'] = 'sa'
            else:
                item['presence'] = 's'

    def gene_experiment(self, gene_id, experiment_names_list):
        """
        Get data to populate javascript splice graph.
        :param gene_id: gene id
        :param experiment_names_list: experiment names
        :return: dictionary
        """

        junc_reads = {}
        ir_reads = {}
        combined_junc_reads = {}
        combined_ir_reads = {}
        all_junctions = list(self.junctions(gene_id, omit_simplified=self.omit_simplified))
        all_ir = list(self.intron_retentions(gene_id, omit_simplified=self.omit_simplified))


        if not experiment_names_list:
            exp_name = "splice graph"
            junc_reads[exp_name] = {}
            ir_reads[exp_name] = {}
            for junc in all_junctions:
                junc_start, junc_end = itemgetter('start', 'end')(junc)
                try:
                    junc_reads[exp_name][junc_start][junc_end] = 0
                except KeyError:
                    junc_reads[exp_name][junc_start] = {junc_end: 0}
            for ir in all_ir:
                ir_start, ir_end = itemgetter('start', 'end')(ir)
                try:
                    ir_reads[exp_name][ir_start][ir_end] = 0
                except KeyError:
                    ir_reads[exp_name][ir_start] = {ir_end: 0}

        else:
            all_junc_reads = self.junction_reads_exp_opt(gene_id)
            all_ir_reads = self.intron_retention_reads_exp_opt(gene_id)
            for experiment_names in experiment_names_list:
                combined_name = next((n for n in experiment_names if ' Combined' in n), '')
                experiment_names = [e for e in experiment_names if e != combined_name]
                combined_junc_reads = {}
                combined_ir_reads = {}

                for name in experiment_names:
                    junc_reads[name] = {}
                    ir_reads[name] = {}
                    if combined_name:
                        junc_reads[combined_name] = {}
                        ir_reads[combined_name] = {}


                for junc in all_junctions:
                    junc_start, junc_end = itemgetter('start', 'end')(junc)

                    for exp_name in experiment_names:
                        try:
                            reads = all_junc_reads[(exp_name, junc_start, junc_end)]
                        except:
                            continue

                        try:
                            junc_reads[exp_name][junc_start][junc_end] = reads
                        except KeyError:
                            junc_reads[exp_name][junc_start] = {junc_end: reads}

                        if combined_name:
                            try:
                                combined_junc_reads[junc_start][junc_end].append(reads)

                            except KeyError:
                                combined_junc_reads[junc_start] = {junc_end: [reads]}

                    if combined_name:

                        try:
                            median_reads = ceil(median(combined_junc_reads[junc_start][junc_end]))
                        except StatisticsError:
                            median_reads = 0
                        except KeyError:
                            median_reads = 0

                        try:
                            junc_reads[combined_name][junc_start][junc_end] = median_reads
                        except KeyError:
                            junc_reads[combined_name][junc_start] = {junc_end: median_reads}

                for ir in all_ir:
                    ir_start, ir_end = itemgetter('start', 'end')(ir)

                    for exp_name in experiment_names:

                        try:
                            reads = all_ir_reads[(exp_name, ir_start, ir_end)]
                        except:
                            continue


                        try:
                            ir_reads[exp_name][ir_start][ir_end] = reads
                            if combined_name:
                                combined_ir_reads[ir_start][ir_end].append(reads)
                        except KeyError:
                            ir_reads[exp_name][ir_start] = {ir_end: reads}
                            if combined_name:
                                combined_ir_reads[ir_start] = {ir_end: [reads]}

                    if combined_name:

                        try:
                            median_reads = ceil(median(combined_ir_reads[ir_start][ir_end]))
                        except StatisticsError:
                            median_reads = 0
                        except KeyError:
                            median_reads = 0

                        try:
                            ir_reads[combined_name][ir_start][ir_end] = median_reads
                        except KeyError:
                            ir_reads[combined_name][ir_start] = {ir_end: median_reads}

        gene_dict = dict(self.view_gene(gene_id))
        gene_dict['exons'] = tuple(dict(e) for e in self.view_exons(gene_id))
        gene_dict['junctions'] = tuple(dict(j) for j in self.view_junctions(gene_id))
        self._add_presence_key(gene_dict['junctions'])
        gene_dict['intron_retention'] = tuple(dict(ir) for ir in self.view_intron_retentions(gene_id))
        self._add_presence_key(gene_dict['intron_retention'])
        gene_dict['junction_reads'] = junc_reads
        gene_dict['intron_retention_reads'] = ir_reads
        gene_dict['genome'] = self.genome
        gene_dict['alt_starts'] = tuple(list(a.values())[0] for a in self.alt_starts(gene_id))
        gene_dict['alt_ends'] = tuple(list(a.values())[0] for a in self.alt_ends(gene_id))

        return gene_dict

    def gene_experiment_combined_only(self, gene_id, experiment_names_list):
        """
        Get data to populate javascript splice graph.
        :param gene_id: gene id
        :param experiment_names_list: experiment names
        :return: dictionary
        """

        junc_reads = {}
        ir_reads = {}
        combined_junc_reads = {}
        combined_ir_reads = {}
        all_junctions = list(self.junctions(gene_id, omit_simplified=self.omit_simplified))
        all_ir = list(self.intron_retentions(gene_id, omit_simplified=self.omit_simplified))


        exp_names = []
        groups = {}
        for experiment_names in experiment_names_list:
            combined_name = next((n for n in experiment_names if ' Combined' in n), '')
            if not combined_name:
                continue
            exp_names.append([combined_name])
            groups[combined_name] = [e for e in experiment_names if e != combined_name]

        all_junc_reads = self.junction_reads_sums(gene_id, groups)
        all_ir_reads = self.intron_retention_reads_sums(gene_id, groups)

        for combined_name in groups.keys():
            junc_reads[combined_name] = {}
            ir_reads[combined_name] = {}


            for junc in all_junctions:
                junc_start, junc_end = itemgetter('start', 'end')(junc)
                median_reads = all_junc_reads[combined_name].get((junc_start, junc_end), 0)
                try:
                    junc_reads[combined_name][junc_start][junc_end] = median_reads
                except KeyError:
                    junc_reads[combined_name][junc_start] = {junc_end: median_reads}


            for ir in all_ir:
                ir_start, ir_end = itemgetter('start', 'end')(ir)
                median_reads = all_ir_reads[combined_name].get((ir_start, ir_end), 0)
                try:
                    ir_reads[combined_name][ir_start][ir_end] = median_reads
                except KeyError:
                    ir_reads[combined_name][ir_start] = {ir_end: median_reads}



        gene_dict = dict(self.view_gene(gene_id))
        gene_dict['exons'] = tuple(dict(e) for e in self.view_exons(gene_id))
        gene_dict['junctions'] = tuple(dict(j) for j in self.view_junctions(gene_id))
        gene_dict['intron_retention'] = tuple(dict(ir) for ir in self.view_intron_retentions(gene_id))
        gene_dict['junction_reads'] = junc_reads
        gene_dict['intron_retention_reads'] = ir_reads
        gene_dict['genome'] = self.genome
        gene_dict['alt_starts'] = tuple(list(a.values())[0] for a in self.alt_starts(gene_id))
        gene_dict['alt_ends'] = tuple(list(a.values())[0] for a in self.alt_ends(gene_id))

        return gene_dict, exp_names

    def gene_transcript_exons(self, gene_id):
        """
        Get all exons for specified gene id
        :param gene_id: gene id
        :return: list of exons
        """
        try:
            query = self.conn.execute('''
                                SELECT gene_id, transcript_id, start, end 
                                FROM transcript_exon 
                                WHERE gene_id=?
                                ORDER BY transcript_id
                                ''', (gene_id,))
        except sqlite3.OperationalError:
            return {}
        transcripts = {}
        cur_id = None
        for row in self._iter_results(query, transcript_exon_fieldnames):
            if row['transcript_id'] != cur_id:
                transcripts[row['transcript_id']] = []
                cur_id = row['transcript_id']
            transcripts[row['transcript_id']].append({'start': row['start'], 'end': row['end'], 'color': 'grey'})
        return transcripts