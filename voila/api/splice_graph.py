import os
import sqlite3
from operator import itemgetter
from pathlib import Path

from voila.constants import EXEC_DIR


class SpliceGraphSQL:
    def __init__(self, filename, delete=False):
        """
        Splice graph class to handle metadata and table information retrieval.
        :param filename: database file name
        :param delete: delete existing database file.
        """

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

        # create tables in freshly delete file.
        if delete is True:
            with (Path(EXEC_DIR) / 'api/model.sql').open('r') as sql:
                self.conn.executescript(sql.read())
                self.conn.commit()

        # verify that file is a sqlite database.
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
        """
        Get genome from data base.
        :return: string
        """
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
        """
        Write genome to databse.
        :param g: genome value
        :return: None
        """
        self.conn.execute('''
                            INSERT 
                            INTO genome (name) 
                            VALUES (?)
                            ''', (g,))

    @property
    def experiment_names(self):
        """
        Get experiment names for database.
        :return: list of strings
        """
        if self._experiment_names is None:
            query = self.conn.execute('''
                                        SELECT name from experiment 
                                        ''')
            fetch = query.fetchall()
            self._experiment_names = tuple(e for e, in fetch)
        return self._experiment_names

    @experiment_names.setter
    def experiment_names(self, names):
        """
        Write experiment names to database
        :param names: list of experiment names
        :return: None
        """
        self.conn.executemany('''
                            INSERT 
                            INTO experiment (name)
                            VALUES (?)
                            ''', tuple((n,) for n in names))

    @property
    def file_version(self):
        """
        Get file version from database.
        :return: string
        """
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
        """
        Write file version from constants to database.
        :param version:
        :return: None
        """
        self.conn.execute('''
                            INSERT 
                            INTO main.file_version (value)
                            VALUES (?)
                            ''', version)

    def _iter_results(self, query, fieldnames):
        """
        Query blocks of 100 entries from database and then convert each line to a dictionary.

        :param query: database query
        :param fieldnames: table fieldnames
        :return: generator of dictionaries
        """
        while True:
            fetch = query.fetchmany(100)
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
alt_starts_fieldnames = ('coordinate',)
alt_ends_fieldnames = ('coordinate',)


class Genes(SpliceGraphSQL):
    def genes(self):
        """
        Get all the genes in the database.
        :return: generator of dictionaries
        """
        query = self.conn.execute('SELECT id, name, strand, chromosome FROM gene')
        return self._iter_results(query, gene_fieldnames)

    def gene(self, gene_id):
        """
        Get gene with gene id supplied.
        :param gene_id: gene id
        :return: dictionary
        """
        query = self.conn.execute('SELECT id, name, strand, chromosome FROM gene WHERE id=?', (gene_id,))
        fetch = query.fetchone()
        if fetch:
            return dict(zip(gene_fieldnames, fetch))


class Exons(SpliceGraphSQL):
    def exons(self, gene_id):
        """
        Get all exons for specified gene id
        :param gene_id: gene id
        :return: list of exons
        """

        query = self.conn.execute('''
                                SELECT gene_id, start, end, annotated_start, annotated_end, annotated 
                                FROM exon 
                                WHERE gene_id=?
                                ''', (gene_id,))
        return self._iter_results(query, exon_fieldnames)


class Junctions(SpliceGraphSQL):
    def junctions(self, gene_id):
        """
        Get list of junctions for specified gene id.
        :param gene_id: gene id
        :return: list of junction dictionaries
        """

        query = self.conn.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM junction 
                                WHERE gene_id=?
                                ''', (gene_id,))
        return self._iter_results(query, junc_fieldnames)

    def junction_reads_exp(self, junction, experiment_names):
        """
        for a junction and a set of experiment names, get a list of reads.
        :param junction: junction dictionary
        :param experiment_names: list of experiment names
        :return: list of reads dictionaries
        """

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
    def intron_retentions(self, gene_id):
        """
        Get all intron retentions for a gene id.
        :param gene_id: gene id
        :return: list of intron retention dictionaries
        """
        query = self.conn.execute('''
                                SELECT gene_id, start, end, has_reads, annotated
                                FROM intron_retention
                                WHERE gene_id=?
                                ''', (gene_id,))

        return self._iter_results(query, ir_fieldnames)

    def intron_retention_reads_exp(self, ir, experiment_names):
        """
        For an intron retention dictionary and a set of experiment names, get a list of intron retention reads.
        :param ir: intron retention dictionary
        :param experiment_names: list of experiment reads.
        :return: list of intron retention reads
        """

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


class AltStarts(SpliceGraphSQL):
    def alt_starts(self, gene_id):
        """
        Get alternate starts for a specific gene.
        :param gene_id: gene id
        :return: list of alt starts dictionary
        """
        query = self.conn.execute('''
                                    SELECT coordinate 
                                    FROM alt_start
                                    WHERE gene_id=?
                                    ''', (gene_id,))
        return self._iter_results(query, alt_starts_fieldnames)


class AltEnds(SpliceGraphSQL):
    def alt_ends(self, gene_id):
        """
        Get alternate ends for specific gene.
        :param gene_id: gene id
        :return: List of alt ends dictionary
        """
        query = self.conn.execute('''
                                    SELECT coordinate 
                                    FROM alt_end
                                    WHERE gene_id=?
                                    ''', (gene_id,))
        return self._iter_results(query, alt_ends_fieldnames)
