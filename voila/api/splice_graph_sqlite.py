import sqlite3
from collections import namedtuple

from voila.api.splice_graph_abstract import SpliceGraphSQLAbstract, GenesAbstract, JunctionsAbstract, \
    IntronRetentionAbstract, ExonsAbstract


class SpliceGraphSQL(SpliceGraphSQLAbstract):
    def __init__(self, filename):
        self.conn = sqlite3.connect(filename)
        self.c = self.conn.cursor()
        self.c.arraysize = 10

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        pass

    @property
    def genome(self):
        pass

    @property
    def experiment_names(self):
        pass

    @property
    def file_version(self):
        pass


Gene = namedtuple('Gene', ('id', 'name', 'strand', 'chromosome'))
Exon = namedtuple('Exon', ('gene_id', 'start', 'end', 'annotated_start', 'annotated_end', 'annotated'))
Junction = namedtuple('Junction', ('gene_id', 'start', 'end', 'has_reads', 'annotated'))
IntronRetention = namedtuple('IntronRetention', ('gene_id', 'start', 'end', 'has_reads', 'annotated'))
JunctionReads = namedtuple('JunctionReads', ('reads', 'experiment_name'))
IntronRetentionReads = namedtuple('IntronRetentionReads', ('reads', 'experiment_name'))


class Genes(GenesAbstract, SpliceGraphSQL):
    def genes(self):
        query = self.c.execute("SELECT id, name, strand, chromosome FROM gene")
        fetch = query.fetchmany()
        while fetch:
            for gene in fetch:
                yield Gene(*gene)
            fetch = query.fetchmany()

    def gene(self, gene_id):
        query = self.c.execute("SELECT id, name, strand, chromosome FROM gene WHERE id=?", (gene_id,))
        fetch = query.fetchone()
        return Gene(*fetch)


class Exons(ExonsAbstract, SpliceGraphSQL):
    def exon(self, gene_id, start, end):
        query = self.c.execute('''
                                SELECT gene_id, start, end, annotated_start, annotated_end, annotated 
                                FROM exon 
                                WHERE gene_id=?
                                AND start=?
                                AND end=?
                                ''', (gene_id, start, end))
        fetch = query.fetchone()
        return Exon(*fetch)

    def exons(self, gene=None):
        if gene:
            query_args = ('''
            SELECT gene_id, start, end, annotated_start, annotated_end, annotated 
            FROM exon 
            WHERE gene_id=?
            ''', (gene.id,))
        else:
            query_args = ("SELECT gene_id, start, end, annotated_start, annotated_end, annotated FROM exon",)
        query = self.c.execute(*query_args)
        fetch = query.fetchmany()
        while fetch:
            for exon in fetch:
                yield Exon(*exon)
            fetch = query.fetchmany()


class Junctions(JunctionsAbstract, SpliceGraphSQL):
    def junction(self, gene_id, start, end):
        query = self.c.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM junction
                                WHERE gene_id=?
                                AND start=?
                                AND end=?
                                ''', (gene_id, start, end))
        fetch = query.fetchone()
        return Junction(*fetch)

    def junctions(self, gene):
        query = self.c.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM junction 
                                WHERE gene_id=?
                                ''', (gene.id,))
        fetch = query.fetchmany()
        while fetch:
            for j in fetch:
                yield Junction(*j)
            fetch = query.fetchmany()

    def junction_reads(self, junction):
        query = self.c.execute('''
                                SELECT reads, experiment_name 
                                FROM junction_reads
                                WHERE junction_gene_id=?
                                AND junction_start=?
                                AND junction_end=?
                                ''', (junction.gene_id, junction.start, junction.end))
        fetch = query.fetchmany()
        while fetch:
            for jr in fetch:
                yield JunctionReads(*jr)
            fetch = query.fetchmany()


class IntronRetentions(IntronRetentionAbstract, SpliceGraphSQL):
    def intron_retention(self, gene_id, start, end):
        query = self.c.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM intron_retention
                                WHERE gene_id=?
                                AND start=?
                                AND end=?
                                ''', (gene_id, start, end))
        fetch = query.fetchone()
        return IntronRetention(*fetch)

    def intron_retentions(self, gene):
        query = self.c.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM intron_retention
                                WHERE gene_id=?
                                ''', (gene.id,))
        fetch = query.fetchmany()
        while fetch:
            for ir in fetch:
                yield IntronRetention(*ir)
            fetch = query.fetchmany()

    def intron_retention_reads(self, intron_retention):
        query = self.c.execute('''
                                SELECT reads, experiment_name 
                                FROM intron_retention_reads
                                WHERE intron_retention_gene_id=?
                                AND intron_retention_start=?
                                AND intron_retention_end=?
                                ''', (intron_retention.gene_id, intron_retention.start, intron_retention.end))
        fetch = query.fetchmany()
        while fetch:
            for iir in fetch:
                yield IntronRetentionReads(*iir)
            fetch = query.fetchmany()
