from voila.config import SplitterConfig, RecombineConfig
import os, sys
import shutil
from voila.voila_log import voila_log
from voila.filter import open_db, copy_table
import h5py
from multiprocessing import Manager, Pool
import time
from voila.api.view_splice_graph import ViewSpliceGraph
import ast

"""
Utils for splitting and recombining a large set of input files for
"""

def split_gene_ids(n):
    with ViewSpliceGraph() as sg:
        all_gene_ids = sg.gene_ids
    for i in range(n):
        yield all_gene_ids[i::n]

def filter_splicegraph(args):
    gene_ids, output_path, q = args
    log = voila_log()
    config = SplitterConfig()
    # sqlite part

    log.debug("Working on %s" % output_path)

    src_conn = open_db(config.splice_graph_file)
    dest_conn = open_db(output_path)

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

    log.debug("Finished %s" % output_path)

    q.put(None)


def filter_voila_file(args):
    voila_file, gene_ids, output_path, q = args
    log = voila_log()

    log.debug("Working on %s" % output_path)

    with h5py.File(voila_file, 'r', libver='latest') as m, h5py.File(output_path, 'w', libver='latest') as m_new:
        # m.lsv_ids()
        main_grp = m_new.create_group('lsvs')
        # new_meta_group = m_new.create_group('metadata')
        m.copy('metadata', m_new)
        for gene_id in gene_ids:
            try:
                m.copy('lsvs/%s' % gene_id, main_grp)
            except KeyError:
                pass

    q.put(None)



def splitter():
    config = SplitterConfig()
    log = voila_log()

    # main strat, extract all gene_ids, split list evenly, run filter stuff
    if not os.path.exists(config.directory):
        os.makedirs(config.directory)

    if os.listdir(config.directory):
        if config.overwrite:
            for _file in os.listdir(config.directory):
                shutil.rmtree(os.path.join(config.directory, _file))
        else:
            voila_log().critical(
                "Files already present in %s; not running classifier "
                "(--overwrite to delete files in this directory and run anyway)" % config.directory)
            sys.exit(1)

    if config.copy_only:

        for i in range(1, config.num_divisions+1):
            part_dir = os.path.join(config.directory, str(i))
            os.makedirs(part_dir, exist_ok=True)

            shutil.copyfile(config.splice_graph_file, os.path.join(part_dir,
                                                                   os.path.basename(config.splice_graph_file)))

            for voila_file in config.voila_files:
                shutil.copyfile(voila_file, os.path.join(part_dir, os.path.basename(voila_file)))

        return

    # voila files part

    manager = Manager()
    q = manager.Queue()

    pool = Pool(config.nproc)


    work_size = 0

    for i, partial_gene_ids in enumerate(split_gene_ids(config.num_divisions), start=1):
        part_dir = os.path.join(config.directory, str(i))
        os.makedirs(part_dir, exist_ok=True)

        output_path = os.path.join(part_dir, os.path.basename(config.splice_graph_file))
        pool.apply_async(filter_splicegraph, [(partial_gene_ids, output_path, q)])
        work_size += 1

        for voila_file in config.voila_files:
            output_path = os.path.join(part_dir, os.path.basename(voila_file))
            pool.apply_async(filter_voila_file, [(voila_file, partial_gene_ids, output_path, q)])
            work_size += 1



    log.info('output %d files' % work_size)

    # monitor loop
    while True:

        size = q.qsize()
        print('Running split [%d/%d]\r' % (size, work_size), end="")
        time.sleep(1)
        if size == work_size:
            break

    #print('                                                  \r', end="")










def recombine():

    log = voila_log()

    tsv_files_to_combine = {}
    config = RecombineConfig()
    dirs = ast.literal_eval(config.directories)
    num_files = 0
    num_unique_files = 0
    for _dir in dirs:
        for root, dirs, files in os.walk(_dir):
            for filename in files:
                if filename.endswith(".tsv"):
                    if not filename in tsv_files_to_combine:
                        num_unique_files += 1
                        tsv_files_to_combine[filename] = []
                    num_files += 1
                    tsv_files_to_combine[filename].append(os.path.join(root, filename))

    log.info("Found %d total TSVs (from %d unique locations) to combine" % (num_files, num_unique_files))
    if not os.path.exists(config.directory):
        os.makedirs(config.directory, exist_ok=True)

    work_size = None
    for i, _tsv in enumerate(tsv_files_to_combine):
        if not work_size:
            work_size = len(tsv_files_to_combine) * len(tsv_files_to_combine[_tsv])
        wrote_headers = False
        with open(os.path.join(config.directory, _tsv), "w") as outfile:
            for j, f in enumerate(tsv_files_to_combine[_tsv]):
                print('Writing %s [%d/%d]                          \r' % (_tsv, i*len(tsv_files_to_combine[_tsv])+j, work_size), end="")
                with open(f, "r") as infile:
                    if wrote_headers is False:
                        wrote_headers = True
                    else:
                        next(infile)
                    for line in infile:
                        outfile.write(line)

                #os.remove(f)

    log.info("Done writing all TSVs!")