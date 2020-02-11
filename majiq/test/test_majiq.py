
import os, sys, tempfile, shutil, time
from subprocess import Popen, PIPE, STDOUT, check_output
from distutils.dir_util import copy_tree
import pandas as pd
from typing import List
import sqlite3
import numpy as np

PRINT_SUBPROC = True

def run_majiq_build_sql(temp_dir):
    """
    Run the build to make the sql file we will compare against
    """

    os.environ['PYTHONPATH'] = ''
    cmd = (
           'majiq',
           'build', 'DB.gff3',
           '-c', 'settings.ini',
           '-j', '4', '-o', 'tmp_build'
           )

    p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
    for line in p.stdout:
        if PRINT_SUBPROC:
            print(line.decode().replace('\n', ''))

def run_majiq_build_majiqs(temp_dir):
    """
    Run the build to make the majiq files we will be comparing
    """

    os.environ['PYTHONPATH'] = ''

    # cmd = (
    #        'majiq',
    #        'build', 'DB.gff3',
    #        '-c', 'settings.ini',
    #        '-j', '4', '-o', 'bam', "--junc-files-only", '--incremental'
    #        )
    #
    # p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
    # for line in p.stdout:
    #     if PRINT_SUBPROC:
    #         print(line.decode().replace('\n', ''))

    # convert paths to absolute

    # replace "bamdirs" with "sjdirs" in config file
    # with open(os.path.join(temp_dir, 'settings.ini')) as f:
    #     newText = f.read().replace('bamdirs', 'sjdirs')
    #
    # with open(os.path.join(temp_dir, 'settings.ini'), "w") as f:
    #     f.write(newText)
    #
    # with open(os.path.join(temp_dir, 'settings.ini')) as f:
    #     print(f.read())

    cmd = (
        'majiq',
        'build', 'DB.gff3',
        '-c', 'settings.ini',
        '-j', '4', '-o', 'tmp_majiqs', '--incremental'
    )

    p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
    for line in p.stdout:
        if PRINT_SUBPROC:
            print(line.decode().replace('\n', ''))


def test_sqlite_equiv(tbl: str, index: List[str], sg_ref: str, sg_test: str) -> None:
    """ Test if specified table is equivalent in sg_ref/sg_test

    Test if specified table is equivalent in sg_ref/sg_test. Ignores/accounts
    for differences in order using unique index

    Parameters
    ----------
    tbl: str
        Name of table in sqlite database we are testing
    index: str
        unique index for tbl
    sg_ref, sg_test: str
        Paths to reference/test splicegraphs (or other database) we are testing
    """
    query = f"select * from {tbl}"
    # load/sort tables
    df_ref: pd.DataFrame = (
        # load table
        pd.read_sql_query(query, sqlite3.connect(sg_ref))
        # set/sort unique index
        .set_index(index, verify_integrity=True).sort_index()
    )
    df_test: pd.DataFrame = (
        # load table
        pd.read_sql_query(query, sqlite3.connect(sg_test))
        # set/sort unique index
        .set_index(index, verify_integrity=True).sort_index()
    )
    # really nice assertion for use with pytest or other Python test frameworks
    pd.testing.assert_frame_equal(df_ref, df_test)
    return

temp_dir = tempfile.mkdtemp()
test_data_dir = os.path.join(os.path.dirname(__file__), 'cases')

def test_workshop1():
    """
    Run a test over the workshop data to determine if the same splicegraph is generated
    Also check that the same majiq files are generated when using the original workshop sjfiles as input
    """

    try:
        shutil.rmtree(temp_dir)
    except:
        pass

    shutil.copytree(os.path.join(test_data_dir, 'workshop/inputs_original/'), temp_dir)


    run_majiq_build_sql(temp_dir)


    sg_ref_path = os.path.join(test_data_dir, 'workshop/outputs_validate/splicegraph.sql')
    sg_test_path = os.path.join(temp_dir, 'tmp_build', 'splicegraph.sql')

    # check splicegraph equalities
    for table_name, table_index in (('alt_end', ["gene_id", "coordinate"]),
                                    ('alt_start', ["gene_id", "coordinate"]),
                                    ('exon', ["gene_id", "start", "end", "annotated_start", "annotated_end", "annotated"]),
                                    ('experiment', ["name"]),
                                    ('file_version', ["id", "value"]),
                                    ('gene', ["id", "name", "strand", "chromosome"]),
                                    ('gene_overlap', ["gene_id_1", "gene_id_2"]),
                                    ('genome', ["id", "name"]),
                                    ('intron_retention', ["gene_id", "start", "end", "has_reads", "annotated", "is_simplified", "is_constitutive"]),
                                    ('intron_retention_reads', ["reads", "experiment_name", "intron_retention_gene_id", "intron_retention_start", "intron_retention_end"]),
                                    ('junction', ["gene_id", "start", "end", "has_reads", "annotated", "is_simplified"]),
                                    ('junction_reads', ["reads", "experiment_name", "junction_gene_id", "junction_start", "junction_end"])):

        test_sqlite_equiv(table_name, table_index, sg_ref_path, sg_test_path)


    shutil.rmtree(temp_dir)

    # move the original data back to the folder
    shutil.copytree(os.path.join(test_data_dir, 'workshop/inputs_original/'), temp_dir)

    # move the partially pre run sj files from the build into the sj dir
    shutil.copytree(os.path.join(test_data_dir, 'workshop/inputs_partial/'), os.path.join(temp_dir, 'sjs'))

    # overwrite the original config file
    os.rename(os.path.join(temp_dir, 'sjs', 'settings.ini'), os.path.join(temp_dir, 'settings.ini'))


    run_majiq_build_majiqs(temp_dir)

    for maj_file in os.listdir(os.path.join(temp_dir, 'tmp_majiqs')):
        if not maj_file.endswith('.majiq'):
            continue

        original_maj_file = os.path.join(os.path.join(test_data_dir, 'workshop/outputs_validate/'), os.path.basename(maj_file))
        testing_maj_file = os.path.join(temp_dir, 'tmp_majiqs', maj_file)

        original_np = np.load(original_maj_file)
        testing_np = np.load(testing_maj_file)

        print(original_maj_file, testing_maj_file)

        for array_name in ('junc_info', 'lsv_types', 'meta', 'coverage',):
            assert np.array_equiv(np.sort(original_np[array_name], axis=0), np.sort(testing_np[array_name], axis=0))


    try:
        shutil.rmtree(temp_dir)
    except:
        pass


"""
Majiq partial-run testing process

majiq build --incremental --junc-files-only  -c /home/pjewell/PycharmProjects/majiq/majiq/test/cases/workshop/inputs_original/settings.ini -o /home/pjewell/PycharmProjects/majiq/majiq/test/cases/workshop/tmp_test/ /home/pjewell/PycharmProjects/majiq/majiq/test/cases/workshop/inputs_original/DB.gff3
majiq build --incremental -c /home/pjewell/PycharmProjects/majiq/majiq/test/cases/workshop/inputs_partial/settings.ini -o /home/pjewell/PycharmProjects/majiq/majiq/test/cases/workshop/outputs/ /home/pjewell/PycharmProjects/majiq/majiq/test/cases/workshop/inputs_original/DB.gff3

"""


test_workshop1()

