import os
import sqlite3
from operator import itemgetter
from pathlib import Path

from voila.api.splice_graph_abstract import SpliceGraphSQLAbstract
from voila.constants import EXEC_DIR


class SpliceGraphSQL(SpliceGraphSQLAbstract):
    def __init__(self, filename, delete=False):
        try:
            filename = filename.decode("utf-8")
        except AttributeError:
            filename = str(filename)

        if delete is True:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

        self.conn = sqlite3.connect(filename)
        self.conn.execute('pragma foreign_keys=ON')

        if delete is True:
            with open(Path(EXEC_DIR) / 'api/model.sql', 'r') as sql:
                self.conn.executescript(sql.read())
                self.conn.commit()

        self.conn.execute('select * from file_version')

        self._genome = None
        self._experiment_names = None
        self._file_version = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.conn.commit()
        self.conn.close()

    @property
    def genome(self):
        if self._genome is None:
            query = self.conn.execute('SELECT name FROM genome')
            genome, = query.fetchone()
            if not genome:
                self._genome = ''
            else:
                self._genome = genome

        return self._genome

    @genome.setter
    def genome(self, g):
        self.conn.execute('''
                            INSERT 
                            INTO genome (name) 
                            VALUES (?)
                            ''', (g,))

    @property
    def experiment_names(self):
        if self._experiment_names is None:
            query = self.conn.execute('''
                                        SELECT name from experiment 
                                        ''')
            fetch = query.fetchall()
            self._experiment_names = tuple(e for e, in fetch)
        return self._experiment_names

    @experiment_names.setter
    def experiment_names(self, names):
        self.conn.executemany('''
                            INSERT 
                            INTO experiment (name)
                            VALUES (?)
                            ''', tuple((n,) for n in names))

    @property
    def file_version(self):
        try:
            if self._file_version is None:
                query = self.conn.execute('''
                                            SELECT value from file_version
                                            ''')
                file_version, = query.fetchone()
                if not file_version:
                    self._file_version = ''
                else:
                    self._file_version = file_version
            return self._file_version
        except TypeError:
            # File version not found in database.
            return -1

    @file_version.setter
    def file_version(self, version):
        self.conn.execute('''
                            INSERT 
                            INTO main.file_version (value)
                            VALUES (?)
                            ''', version)

    def _iter_results(self, query, fieldnames):
        while True:
            fetch = query.fetchmany(1)
            if not fetch:
                break
            for x in fetch:
                yield dict(zip(fieldnames, x))


gene_fieldnames = ('id', 'name', 'strand', 'chromosome')
junc_fieldnames = ('gene_id', 'start', 'end', 'has_reads', 'annotated')
junc_reads_fieldnames = ('reads', 'experiment_name')
exon_fieldnames = ('gene_id', 'start', 'end', 'annotated_start', 'annotated_end', 'annotated')
ir_fieldnames = ('gene_id', 'start', 'end', 'has_reads', 'annotated')
ir_reads_fieldnames = ('reads', 'experiment_name')


class Genes(SpliceGraphSQL):
    def genes(self):
        query = self.conn.execute('SELECT id, name, strand, chromosome FROM gene')
        return self._iter_results(query, gene_fieldnames)

    def gene(self, gene_id):
        query = self.conn.execute('SELECT id, name, strand, chromosome FROM gene WHERE id=?', (gene_id,))
        fetch = query.fetchone()
        if fetch:
            return dict(zip(gene_fieldnames, fetch))


class Exons(SpliceGraphSQL):
    def exons(self, gene):
        if isinstance(gene, str):
            gene_id = gene
        else:
            gene_id = gene['id']

        query = self.conn.execute('''
                                SELECT gene_id, start, end, annotated_start, annotated_end, annotated 
                                FROM exon 
                                WHERE gene_id=?
                                ''', (gene_id,))
        return self._iter_results(query, exon_fieldnames)


class Junctions(SpliceGraphSQL):
    def junctions(self, gene):
        if isinstance(gene, str):
            gene_id = gene
        else:
            gene_id = gene['id']

        query = self.conn.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM junction 
                                WHERE gene_id=?
                                ''', (gene_id,))
        return self._iter_results(query, junc_fieldnames)

    def junction_reads_exp(self, junction, experiment_names):
        query = self.conn.execute('''
                                SELECT reads, experiment_name 
                                FROM junction_reads
                                WHERE junction_start=?
                                AND junction_end=?
                                AND junction_gene_id=?
                                AND experiment_name IN ({})
                                '''.format(','.join(["'{}'".format(x) for x in experiment_names])),
                                  itemgetter('start', 'end', 'gene_id')(junction))
        return self._iter_results(query, junc_reads_fieldnames)


class IntronRetentions(SpliceGraphSQL):
    def intron_retentions(self, gene):
        query = self.conn.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM intron_retention
                                WHERE gene_id=?
                                ''', (gene['id'],))
        return self._iter_results(query, ir_fieldnames)

    def intron_retention_reads_exp(self, ir, experiment_names):
        query = self.conn.execute('''
                                SELECT reads, experiment_name 
                                FROM intron_retention_reads
                                WHERE intron_retention_start=?
                                AND intron_retention_end=?
                                AND intron_retention_gene_id=?
                                AND experiment_name IN ({})
                                '''.format(','.join(["'{}'".format(x) for x in experiment_names])),
                                  (ir['start'], ir['end'], ir['gene_id']))
        return self._iter_results(query, ir_reads_fieldnames)
