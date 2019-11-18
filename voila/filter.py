from voila.config import FilterConfig
from voila import constants
from voila.voila_log import voila_log
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_hdf5 import MatrixHdf5
from voila.api.view_matrix import ViewHeterogen, ViewPsi, ViewDeltaPsi
from voila.constants import *
import h5py
import sqlite3
import os, sys
import numpy as np
from voila.vlsv import get_expected_psi, matrix_area
from itertools import combinations
from multiprocessing import Manager, Pool
import time
from voila.api.matrix_utils import generate_means, unpack_bins, generate_high_probability_non_changing

class Filter:
    def __init__(self):

        config = FilterConfig()
        analysis_type = config.analysis_type

        voila_log().info(analysis_type + ' FILTER')

        run_filter()

def lsv_ids2gene_ids(lsv_ids):
    gene_ids = set()
    for lsv_id in lsv_ids:
        parts = lsv_id.split(':s:')
        if not len(parts) == 2:
            parts = lsv_id.split(':t:')
        gene_ids.add(parts[0])
    return gene_ids

def read_ids_file(file_path):
    ids = set()
    with open(file_path, 'r') as f:
        for _id in f:
            if not _id.isspace():
                ids.add(_id.replace('\n', ''))
    return ids

def open_db(db_file_path):
    conn = sqlite3.connect(db_file_path)
    # conn.row_factory = sqlite3.Row
    return conn


def divide_chunks(l, n):
    # read from list and return chunks of length n
    for i in range(0, len(l), n):
        yield l[i:i + n]

def copy_table(table, src, dest, gene_ids=None, gene_ids_colname=None):
    if gene_ids and gene_ids_colname:
        scs = []
        for gene_ids_chunk in divide_chunks(gene_ids, 999):
            format_strings = ','.join(['?'] * len(gene_ids_chunk))
            scs.append(src.execute('SELECT * FROM %s WHERE %s in (%s)' % (table, gene_ids_colname, format_strings),
                             tuple(gene_ids_chunk)))
    else:
        scs = [src.execute('SELECT * FROM %s' % table)]
    dc = dest.cursor()

    cols = None
    for i, sc in enumerate(scs):
        if i == 0:
            cols = [description[0] for description in sc.description]
            dc.execute("CREATE TABLE %s (%s)" % (table, ','.join(cols)))
        ins = 'INSERT INTO %s (%s) VALUES (%s)' % (table, ','.join(cols), ','.join(['?'] * len(cols)))
        dc.executemany(ins, sc.fetchall())

    dest.commit()

def filter_splicegraph(args):
    gene_ids, q = args
    voila_log().info("Filtering Splicegraph")
    config = FilterConfig()
    # sqlite part



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

    voila_log().info("Finished Filtering Splicegraph")
    q.put(None)

def filter_voila_file(args):
    voila_file, gene_ids, lsv_ids, q = args
    voila_log().info("Started Filtering %s" % voila_file)
    config = FilterConfig()

    new_voila_file = os.path.join(config.directory, os.path.basename(voila_file))
    if os.path.exists(new_voila_file):
        if config.overwrite:
            os.remove(new_voila_file)
        else:
            voila_log().warning(
                "%s already exists, skipping writing this file. (--overwrite to bypass)" % new_voila_file)
            return

    with MatrixHdf5(voila_file) as m:
        analysis_type = m.analysis_type

    included_lsv_ids = set()
    excluded_lsv_ids = set()
    # for PSI filter, any input file can be used
    # for dPSI filter, only matters for dPSI and HET files

    # if no relevant filters provided, basically skip it
    if not config.changing and not config.non_changing:
            #and not config.decomplexify_psi_threshold > 0 and not config.decomplexify_deltapsi_threshold > 0:
        included_lsv_ids = lsv_ids

    if config.changing:

        if analysis_type == ANALYSIS_PSI:
            voila_log().debug("Skipping changing filter for %s because no delta-psi available" % voila_file)

        elif analysis_type == ANALYSIS_DELTAPSI:
            with ViewDeltaPsi(voila_file) as m:
                if not lsv_ids:
                    lsv_ids = list(m.lsv_ids())

                for lsv in m.lsvs():
                    if lsv.lsv_id in lsv_ids:
                        included_lsv_ids.add(lsv.lsv_id)

                        # calc confidence thresh
                        bins = lsv.bins
                        conf_change = (matrix_area(b, config.changing_threshold) for b in bins)

                        # print(list(v for v in excl_incl))

                        if not any(cc >= config.probability_changing_threshold for cc in conf_change):
                            excluded_lsv_ids.add(lsv.lsv_id)

        elif analysis_type == ANALYSIS_HETEROGEN:
            with ViewHeterogen(voila_file) as m:
                if not lsv_ids:
                    lsv_ids = list(m.lsv_ids())

                for lsv in m.lsvs():
                    if lsv.lsv_id in lsv_ids:
                        included_lsv_ids.add(lsv.lsv_id)

                        if not lsv.lsv_id in excluded_lsv_ids:

                            conf_change = (matrix_area(b, config.changing_threshold) for b in generate_means(np.array(lsv.mean_psi).transpose((1, 0, 2))))
                            if not any(cc >= config.probability_changing_threshold for cc in conf_change):
                                excluded_lsv_ids.add(lsv.lsv_id)

    elif config.non_changing:
        if analysis_type == ANALYSIS_PSI:
            voila_log().debug("Skipping changing filter for %s because no delta-psi available" % voila_file)

        elif analysis_type == ANALYSIS_DELTAPSI:
            with ViewDeltaPsi(voila_file) as m:
                if not lsv_ids:
                    lsv_ids = list(m.lsv_ids())

                for lsv in m.lsvs():
                    if lsv.lsv_id in lsv_ids:
                        included_lsv_ids.add(lsv.lsv_id)

                        # calc confidence thresh


                        # print(list(v for v in excl_incl))

                        conf_nonchange = max(generate_high_probability_non_changing(lsv.intron_retention,
                                                                           lsv.matrix_hdf5.prior,
                                                                           config.non_changing_threshold,
                                                                           lsv.bins))

                        #if not all(cc >= config.probability_non_changing_threshold for cc in conf_nonchange):
                        if not conf_nonchange >= config.probability_non_changing_threshold:
                            excluded_lsv_ids.add(lsv.lsv_id)

        elif analysis_type == ANALYSIS_HETEROGEN:
            with ViewHeterogen(voila_file) as m:
                if not lsv_ids:
                    lsv_ids = list(m.lsv_ids())

                for lsv in m.lsvs():
                    if lsv.lsv_id in lsv_ids:
                        included_lsv_ids.add(lsv.lsv_id)

                        if not lsv.lsv_id in excluded_lsv_ids:



                            # conf_change = (matrix_area(b, config.changing_threshold) for b in
                            #                generate_means(np.array(lsv.mean_psi).transpose((1, 0, 2))))
                            # if not any(cc >= config.probability_changing_threshold for cc in conf_change):
                            if False:
                                excluded_lsv_ids.add(lsv.lsv_id)


    # if config.decomplexify_psi_threshold > 0 or config.decomplexify_deltapsi_threshold > 0:
    #     # standard psi / dpsi based filters
    #
    #     if analysis_type == ANALYSIS_PSI:
    #         with ViewPsi(voila_file) as m:
    #             if not lsv_ids:
    #                 lsv_ids = list(m.lsv_ids())
    #
    #             for lsv in m.lsvs():
    #                 if lsv.lsv_id in lsv_ids:
    #                     included_lsv_ids.add(lsv.lsv_id)
    #                     if config.decomplexify_psi_threshold > 0:
    #                         if any(v < config.decomplexify_psi_threshold for v in lsv.means):
    #                             excluded_lsv_ids.add(lsv.lsv_id)
    #                             continue
    #
    #     elif analysis_type == ANALYSIS_DELTAPSI:
    #         with ViewDeltaPsi(voila_file) as m:
    #             if not lsv_ids:
    #                 lsv_ids = list(m.lsv_ids())
    #
    #             for lsv in m.lsvs():
    #                 if lsv.lsv_id in lsv_ids:
    #                     included_lsv_ids.add(lsv.lsv_id)
    #                     if config.decomplexify_psi_threshold > 0:
    #                         for means in lsv.group_means:
    #                             # means: (name, means,)
    #                             if any(v < config.decomplexify_psi_threshold for v in means[1]):
    #                                 excluded_lsv_ids.add(lsv.lsv_id)
    #                                 break
    #
    #                     if config.decomplexify_deltapsi_threshold > 0:
    #                         excl_incl = lsv.excl_incl
    #                         # print(list(v for v in excl_incl))
    #                         if any((v[0] > 0 and v[0] < config.decomplexify_deltapsi_threshold) or \
    #                                (v[1] > 0 and v[1] < config.decomplexify_deltapsi_threshold) for v in excl_incl):
    #                             excluded_lsv_ids.add(lsv.lsv_id)
    #
    #
    #
    #     elif analysis_type == ANALYSIS_HETEROGEN:
    #         with ViewHeterogen(voila_file) as m:
    #             if not lsv_ids:
    #                 lsv_ids = list(m.lsv_ids())
    #
    #             for lsv in m.lsvs():
    #                 if lsv.lsv_id in lsv_ids:
    #                     included_lsv_ids.add(lsv.lsv_id)
    #                     if config.decomplexify_psi_threshold > 0:
    #                         if not lsv.lsv_id in excluded_lsv_ids:
    #                             for mean in np.array(lsv.mean_psi).transpose((1, 0, 2)):
    #                                 if any(get_expected_psi(v) < config.decomplexify_psi_threshold for v in mean):
    #                                     excluded_lsv_ids.add(lsv.lsv_id)
    #                                     break
    #                     if config.decomplexify_deltapsi_threshold > 0:
    #                         if not lsv.lsv_id in excluded_lsv_ids:
    #                             for psis_g1, psis_g2 in combinations(np.array(lsv.mean_psi).transpose((1, 0, 2)), 2):
    #                                 for psi_g1, psi_g2 in zip(psis_g1, psis_g2):
    #                                     if abs(get_expected_psi(psi_g1) - get_expected_psi(
    #                                             psi_g2)) < config.decomplexify_deltapsi_threshold:
    #                                         excluded_lsv_ids.add(lsv.lsv_id)
    #                                         break
    #                                 else:
    #                                     continue
    #                                 break


    if not gene_ids:
        _gene_ids = lsv_ids2gene_ids(included_lsv_ids)
    else:
        _gene_ids = gene_ids

    if config.changing or config.non_changing:
        voila_log().info("Filtered to %d LSVs for %s, writing output" % (len(included_lsv_ids) - len(excluded_lsv_ids), voila_file))
        if not included_lsv_ids:
            voila_log().critical("For voila file %s, no lsvs matched the specified filters!" % voila_file)
            q.put([])
            return []

    final_lsv_ids = []
    with h5py.File(voila_file, 'r', libver='latest') as m, h5py.File(new_voila_file, 'w', libver='latest') as m_new:
        # m.lsv_ids()
        main_grp = m_new.create_group('lsvs')
        # new_meta_group = m_new.create_group('metadata')
        m.copy('metadata', m_new)
        for gene_id in _gene_ids:

            # new_gene_group = m_new.create_group('lsvs/%s' % gene_id)
            if not lsv_ids and not config.changing and not config.non_changing:
                #and not config.decomplexify_psi_threshold > 0 and not config.decomplexify_deltapsi_threshold > 0:\
                # this part only makes sense if copying ALL lsv ids...
                m.copy('lsvs/%s' % gene_id, main_grp)
            else:
                lsv_grp = m_new.create_group('lsvs/%s' % gene_id)

                for lsv_id in included_lsv_ids:
                    if lsv_id.startswith(gene_id) and not lsv_id in excluded_lsv_ids:
                        final_lsv_ids.append(lsv_id)
                        m.copy('lsvs/%s/%s' % (gene_id, lsv_id), lsv_grp)

    voila_log().info("Finished writing %s" % new_voila_file)


    q.put(final_lsv_ids, True)
    return final_lsv_ids

def run_filter():

    config = FilterConfig()

    num_primary_filters = sum(bool(x) for x in (config.gene_ids, config.gene_ids_file, config.lsv_ids, config.lsv_ids_file))
    if num_primary_filters == 0 and not config.changing and not config.non_changing:
        voila_log().critical(
            "In order to filter, you must specify --gene-ids, --gene-ids-file, --lsv-ids, --lsv-ids-file")
        sys.exit(2)
    if num_primary_filters > 1:
        voila_log().critical(
            "Please specify only one of --gene-ids, --gene-ids-file, --lsv-ids, --lsv-ids-file")
        sys.exit(2)

    if config.lsv_ids:
        lsv_ids = set(config.lsv_ids)
        gene_ids = lsv_ids2gene_ids(lsv_ids)
    elif config.lsv_ids_file:
        lsv_ids = set(read_ids_file(config.lsv_ids_file))
        gene_ids = lsv_ids2gene_ids(lsv_ids)
    elif config.gene_ids:
        lsv_ids = set()
        gene_ids = config.gene_ids
    elif config.gene_ids_file:
        lsv_ids = set()
        gene_ids = read_ids_file(config.gene_ids_file)
    else:
        # only secondary filters
        gene_ids = []
        lsv_ids = set()

    if not os.path.exists(config.directory):
        os.makedirs(config.directory)

    if gene_ids or lsv_ids:
        voila_log().info("Filtering input to %d gene(s) %s" % (len(gene_ids), "; %d LSV(s)" % len(lsv_ids) if lsv_ids else ''))
    voila_log().info("One splicegraph and %d voila files" % len(config.voila_files))
    voila_log().info("Writing output files to %s" % os.path.abspath(config.directory))

    # voila files part

    manager = Manager()
    q = manager.Queue()

    p = Pool(config.nproc)

    # voila_index = p.map(self._heterogen_pool_add_index, zip(lsv_ids, range(work_size), repeat(work_size)))

    final_lsv_ids = set()

    if not config.splice_graph_only:
        voila_files_pool = p.map_async(filter_voila_file,
                                   ((voila_file, gene_ids, lsv_ids, q) for voila_file in config.voila_files), )


        # monitor loop
        while True:

            if voila_files_pool.ready():
                break
            else:
                #size = q.qsize()
                #print('Processing Genes and Modules [%d/%d]\r' % (size, work_size), end="")
                time.sleep(1)

        #print('                                                  \r', end="")
        for voila_file_res in voila_files_pool.get(True):
            for lsv_id in voila_file_res:
                final_lsv_ids.add(lsv_id)

        gene_ids = lsv_ids2gene_ids(final_lsv_ids)

    if not config.voila_files_only:
        splicegraph_pool = p.map_async(filter_splicegraph, ((list(gene_ids), q),))

        # monitor loop
        while True:

            if splicegraph_pool.ready():
                break
            else:
                # size = q.qsize()
                # print('Processing Genes and Modules [%d/%d]\r' % (size, work_size), end="")
                time.sleep(1)

        # print('                                                  \r', end="")
        splicegraph_pool.get()


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