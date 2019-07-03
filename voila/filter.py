from voila.config import FilterConfig
from voila import constants
from voila.voila_log import voila_log
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_hdf5 import MatrixHdf5
import h5py
import sqlite3
import os, sys

class Filter:
    def __init__(self):

        config = FilterConfig()
        analysis_type = config.analysis_type

        voila_log().info(analysis_type + ' FILTER')

        run_filter()


def filter_splicegraph(gene_ids):
    config = FilterConfig()
    # sqlite part

    def open_db(db_file_path):
        conn = sqlite3.connect(db_file_path)
        # conn.row_factory = sqlite3.Row
        return conn

    def copy_table(table, src, dest, gene_ids=None, gene_ids_colname=None):
        if gene_ids and gene_ids_colname:
            format_strings = ','.join(['?'] * len(gene_ids))
            sc = src.execute('SELECT * FROM %s WHERE %s in (%s)' % (table, gene_ids_colname, format_strings),
                             tuple(gene_ids))
        else:
            sc = src.execute('SELECT * FROM %s' % table)
        dc = dest.cursor()

        cols = [description[0] for description in sc.description]
        dc.execute("CREATE TABLE %s (%s)" % (table, ','.join(cols)))
        for row in sc.fetchall():
            ins = 'INSERT INTO %s (%s) VALUES (%s)' % (table, ','.join(cols),
                                                       ','.join(['?'] * len(cols)))

            dc.execute(ins, row)

        dest.commit()

    new_splice_graph_file = os.path.join(config.directory, os.path.basename(config.splice_graph_file))

    if os.path.exists(new_splice_graph_file):
        if config.overwrite:
            os.remove(new_splice_graph_file)
        else:
            voila_log().warning(
                "%s already exists, skipping writing this file. (--overwrite to bypass)" % new_splice_graph_file)
            return

    src_conn = open_db(config.splice_graph_file)
    dest_conn = open_db(new_splice_graph_file)

    # these tables we copy without modifying for now
    copy_table('genome', src_conn, dest_conn)
    copy_table('experiment', src_conn, dest_conn)
    copy_table('file_version', src_conn, dest_conn)
    copy_table('gene_overlap', src_conn, dest_conn)

    # filtered tables
    copy_table('alt_end', src_conn, dest_conn, gene_ids, 'gene_id')
    copy_table('alt_start', src_conn, dest_conn, gene_ids, 'gene_id')
    copy_table('exon', src_conn, dest_conn, gene_ids, 'gene_id')
    copy_table('gene', src_conn, dest_conn, gene_ids, 'id')
    copy_table('intron_retention', src_conn, dest_conn, gene_ids, 'gene_id')
    copy_table('intron_retention_reads', src_conn, dest_conn, gene_ids, 'intron_retention_gene_id')
    copy_table('junction', src_conn, dest_conn, gene_ids, 'gene_id')
    copy_table('junction_reads', src_conn, dest_conn, gene_ids, 'junction_gene_id')


def run_filter():

    config = FilterConfig()

    if not config.gene_ids and not config.gene_ids_file:
        voila_log().critical("In order to filter, you must specify either --gene-ids or --gene-ids-file")
        sys.exit(2)
    elif config.gene_ids and config.gene_ids_file:
        voila_log().critical("In order to filter, you must specify either --gene-ids or --gene-ids-file (not both)")
        sys.exit(2)
    elif config.gene_ids:
        gene_ids = config.gene_ids
    elif config.gene_ids_file:
        gene_ids = []
        with open(config.gene_ids_file, 'r') as f:
            for gene_id in f:
                if not gene_id.isspace():
                    gene_ids.append(gene_id.replace('\n', ''))
    else:
        voila_log().error("Should not get here")
        gene_ids = []



    if not os.path.exists(config.directory):
        os.makedirs(config.directory)

    voila_log().info("Filtering input to %d gene(s)" % len(gene_ids))
    voila_log().info("One splicegraph and %d voila files" % len(config.voila_files))
    voila_log().info("Writing output files to %s" % os.path.abspath(config.directory))

    # voila files part
    for voila_file in config.voila_files:
        new_voila_file = os.path.join(config.directory, os.path.basename(voila_file))
        if os.path.exists(new_voila_file):
            if config.overwrite:
                os.remove(new_voila_file)
            else:
                voila_log().warning("%s already exists, skipping writing this file. (--overwrite to bypass)" % new_voila_file)
                continue

        with h5py.File(voila_file, 'r', libver='latest') as m, h5py.File(new_voila_file, 'w', libver='latest') as m_new:
            #m.lsv_ids()
            lsv_grp = m_new.create_group('lsvs')
            #new_meta_group = m_new.create_group('metadata')
            m.copy('metadata', m_new)
            for gene_id in gene_ids:

                #new_gene_group = m_new.create_group('lsvs/%s' % gene_id)
                m.copy('lsvs/%s' % gene_id, lsv_grp)



    filter_splicegraph(gene_ids)










    voila_log().info("Filtering Complete")


# import multiprocessing
#
#
# def Writer(dest_filename, some_queue, some_stop_token):
#     with open(dest_filename, 'w') as dest_file:
#         while True:
#             line = some_queue.get()
#             if line == some_stop_token:
#                 return
#             dest_file.write(line)
#
#
# def the_job(some_queue):
#     for item in something:
#         result = process(item)
#         some_queue.put(result)
#
#
# if __name__ == "__main__":
#     queue = multiprocessing.Queue()
#
#     STOP_TOKEN = "STOP!!!"
#
#     writer_process = multiprocessing.Process(target=Writer, args=("output.txt", queue, STOP_TOKEN))
#     writer_process.start()
#
#     # Dispatch all the jobs
#
#     # Make sure the jobs are finished
#
#     queue.put(STOP_TOKEN)
#     writer_process.join()
#     # There, your file was written.