import sqlite3
from collections import namedtuple

from voila.api.splice_graph_abstract import SpliceGraphSQLAbstract


class SpliceGraphSQL(SpliceGraphSQLAbstract):
    def __init__(self, filename):
        self.conn = sqlite3.connect(filename)
        self._genome = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.conn.close()

    @property
    def genome(self):
        if not self._genome:
            query = self.conn.execute('''
                                    SELECT name FROM genome
                                    ''')
            genome, = query.fetchone()
            self._genome = genome
        return self._genome

    @property
    def experiment_names(self):
        pass

    @property
    def file_version(self):
        pass

    def _iter_results(self, query, sg_type):
        while True:
            fetch = query.fetchmany(10000)
            if not fetch:
                break
            yield from (sg_type(*x) for x in fetch)


Gene = namedtuple('Gene', ('id', 'name', 'strand', 'chromosome'))
Exon = namedtuple('Exon', ('gene_id', 'start', 'end', 'annotated_start', 'annotated_end', 'annotated'))
Junction = namedtuple('Junction', ('gene_id', 'start', 'end', 'has_reads', 'annotated'))
IntronRetention = namedtuple('IntronRetention', ('gene_id', 'start', 'end', 'has_reads', 'annotated'))
JunctionReads = namedtuple('JunctionReads', ('reads', 'experiment_name'))
IntronRetentionReads = namedtuple('IntronRetentionReads', ('reads', 'experiment_name'))


class Genes(SpliceGraphSQL):
    def genes(self):
        query = self.conn.execute('SELECT id, name, strand, chromosome FROM gene')
        return self._iter_results(query, Gene)

    def gene(self, gene_id):
        query = self.conn.execute('SELECT id, name, strand, chromosome FROM gene WHERE id=?', (gene_id,))
        fetch = query.fetchone()
        return Gene(*fetch)


class Exons(SpliceGraphSQL):
    def exons(self, gene):
        query = self.conn.execute('''
                                SELECT gene_id, start, end, annotated_start, annotated_end, annotated 
                                FROM exon 
                                WHERE gene_id=?
                                ''', (gene.id,))
        return self._iter_results(query, Exon)


class Junctions(SpliceGraphSQL):
    def junctions(self, gene):
        query = self.conn.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM junction 
                                WHERE gene_id=?
                                ''', (gene.id,))
        return self._iter_results(query, Junction)

    def junction_reads(self, junction):
        query = self.conn.execute('''
                                SELECT reads, experiment_name 
                                FROM junction_reads
                                WHERE junction_gene_id=?
                                AND junction_start=?
                                AND junction_end=?
                                ''', (junction.gene_id, junction.start, junction.end))
        return self._iter_results(query, JunctionReads)


class IntronRetentions(SpliceGraphSQL):
    def intron_retentions(self, gene):
        query = self.conn.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM intron_retention
                                WHERE gene_id=?
                                ''', (gene.id,))
        return self._iter_results(query, IntronRetention)

    def intron_retention_reads(self, intron_retention):
        query = self.conn.execute('''
                                SELECT reads, experiment_name 
                                FROM intron_retention_reads
                                WHERE intron_retention_gene_id=?
                                AND intron_retention_start=?
                                AND intron_retention_end=?
                                ''', (intron_retention.gene_id, intron_retention.start, intron_retention.end))
        return self._iter_results(query, IntronRetentionReads)
