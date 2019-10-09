from bisect import bisect_left, bisect_right
from itertools import combinations
from pathlib import Path

from voila import constants
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_utils import generate_means

from voila.classifier.as_types import Graph
from voila.classifier.tsv_writer import TsvWriter


from subprocess import Popen, PIPE, STDOUT
import os, shutil
import csv



# this should vary depending on the group of tests to run
# changed relative import path here:
from tests_expected_t_cells_2 import *
#from tests_expected_t_cells_3 import *




#out_dir = '/Users/calebradens/Documents/majiq_dev/classifier_dev/classified'
out_dir = '/home/paul/PycharmProjects/majiq/test_cases/classifier/caleb1/testout'


def run_voila_classify(gene_ids, enabled_outputs='all', additional_args=[]):
    #os.environ['PYTHONPATH'] = '/Users/calebradens/PycharmProjects/classifier_caleb_dev/'
    os.environ['PYTHONPATH'] = '/home/paul/PycharmProjects/majiq'
    cmd = ['python3',
           #'/Users/calebradens/PycharmProjects/classifier_caleb_dev/voila/run_voila.py',
           '/home/paul/PycharmProjects/majiq/voila/run_voila.py',
           'classify', psi_file, sg_file, '-d', out_dir,
           '--enabled-outputs', enabled_outputs, '--overwrite',
           '--decomplexify-psi-threshold', '0.0']
    for arg in additional_args:
        cmd.append(arg)

    cmd.append('--gene-ids')
    cmd += gene_ids
    #print(' '.join(cmd))
    p = Popen(cmd, stdout=PIPE, stderr=STDOUT)
    output = ''
    error = False
    for line in p.stdout:
        output += line.decode()
        if 'Traceback' in output:
            error = True
        if True:
            print(line.decode().replace('\n', ''))
    if error:
        print(output)
        assert False
    os.environ['PYTHONPATH'] = ''


expected_headers = ['module_id', 'gene_id', 'gene_name', "Chr","Strand", 'lsv_ids', 'cassette_exon', 'tandem_cassette', 'alt3ss', 'alt5ss', 'p_alt3ss',
                    'p_alt5ss', 'alt3and5ss', 'mutually_exclusive', 'alternative_intron', 'ale', 'afe', 'p_ale', 'p_afe', 'orphan_junction',
                    'multi_exon_spanning',   'exitron', 'complex', 'number-of-events']
expected_headers_constitutive = ['module_id', 'gene_id', 'gene_name', "Chr","Strand", 'lsv_ids', 'cassette_exon', 'tandem_cassette', 'alt3ss', 'alt5ss', 'p_alt3ss',
                    'p_alt5ss', 'alt3and5ss', 'mutually_exclusive', 'alternative_intron', 'ale', 'afe', 'p_ale', 'p_afe', 'orphan_junction',
                                 'constitutive_junction', 'constitutive_intron',
                    'multi_exon_spanning',   'exitron', 'complex', 'number-of-events']
expected_headers_mpe= ['Module ID', 'Gene ID', 'Gene Name', "Chr","Strand", 'LSV ID(s)',
                       "Collapsed Event Name","Type","Edge of the Module",
                       "Reference Exon Coord","Reference Exon De Novo","Reference Exon Exitrons","Reference Exon Constant Region",
                       "Reference Exon Trimmed","Constitutive Direction","Constitutive Regions",
                       "Constitutive De Novo","Constitutive Exon or Intron"]


def verify_tsvs(gene_id):

    # with open(os.path.join(out_dir, 'cassette.tsv'), 'r', newline='') as csvfile:
    #     reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')
    #     headers = next(reader, None)
    #     events = []
    #     for line in reader:
    #         if line[1] == gene_id:
    #             events.append(line)
    #
    #     if gene_id in expected_cassette_exons:
    #         for i, mod in enumerate(expected_cassette_exons[gene_id]):
    #             print(mod)
    #
    #             for k, v in mod.items():
    #                 try:
    #                     assert v == events[i][headers.index(k)]
    #                 except:
    #                     print("expt: %s found: %s (%d, %s, %s)" % (v, events[i][headers.index(k)], i+1,
    #                                                                events[i][headers.index(k)], gene_id))
    #                     raise
    #
    # with open(os.path.join(out_dir, 'alt5prime.tsv'), 'r', newline='') as csvfile:
    #     reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')
    #     headers = next(reader, None)
    #     events = []
    #     for line in reader:
    #         if line[1] == gene_id:
    #             events.append(line)
    #
    #     if gene_id in expected_alt5ss:
    #         for i, mod in enumerate(expected_alt5ss[gene_id]):
    #             print(mod)
    #
    #             for k, v in mod.items():
    #                 try:
    #                     assert v == events[i][headers.index(k)]
    #                 except:
    #                     print("expt: %s found: %s (%d, %s, %s)" % (v, events[i][headers.index(k)], i+1,
    #                                                                events[i][headers.index(k)], gene_id))
    #                     raise
    #
    # with open(os.path.join(out_dir, 'alt3prime.tsv'), 'r', newline='') as csvfile:
    #     reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')
    #     headers = next(reader, None)
    #     events = []
    #     for line in reader:
    #         if line[1] == gene_id:
    #             events.append(line)
    #
    #     if gene_id in expected_alt3ss:
    #         for i, mod in enumerate(expected_alt3ss[gene_id]):
    #             print(mod)
    #
    #             for k, v in mod.items():
    #                 try:
    #                     assert v == events[i][headers.index(k)]
    #                 except:
    #                     print("expt: %s found: %s (%d, %s, %s)" % (v, events[i][headers.index(k)], i+1,
    #                                                                events[i][headers.index(k)], gene_id))
    #
    # with open(os.path.join(out_dir, 'alternative_intron.tsv'), 'r', newline='') as csvfile:
    #     reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')
    #     headers = next(reader, None)
    #     events = []
    #     for line in reader:
    #         if line[1] == gene_id:
    #             events.append(line)
    #
    #     if gene_id in expected_alternative_intron:
    #         for i, mod in enumerate(expected_alternative_intron[gene_id]):
    #             print(mod)
    #
    #             for k, v in mod.items():
    #                 try:
    #                     assert v == events[i][headers.index(k)]
    #                 except:
    #                     print("expt: %s found: %s (%d, %s, %s)" % (v, events[i][headers.index(k)], i+1,
    #                                                                events[i][headers.index(k)], gene_id))
    #                     raise

    with open(os.path.join(out_dir, 'summary.tsv'), 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')

        headers = next(reader, None)
        modules = []
        for line in reader:
            if line[1] == gene_id:
                modules.append(line)

        if gene_id in expected_modules:
            if expected_modules[gene_id]:
                try:
                    assert len(modules) == len(expected_modules[gene_id])
                except:
                    if gene_id == "gene:ENSG00000082074":
                        for mod in modules:
                            print(mod)
                    print("expt: %d found: %d (%s)" % (len(expected_modules[gene_id]), len(modules), gene_id))
                    raise


                for i, mod in enumerate(expected_modules[gene_id]):
                    print(modules[i])
                    print(mod)

                    for k, v in mod.items():
                        try:
                            if k == 'lsv_ids':
                                assert all(x in v.split(';') for x in modules[i][expected_headers.index(k)].split(';'))
                            else:

                                assert v == modules[i][expected_headers.index(k)]
                        except:
                            print("expt: %s found: %s (%d, %s, %s)" % (v, modules[i][expected_headers.index(k)], i+1,
                                                                       headers[expected_headers.index(k)], gene_id))
                            raise


def verify_constitutive(gene_id):
    with open(os.path.join(out_dir, 'summary.tsv'), 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')

        headers = next(reader, None)
        modules = []
        for line in reader:
            if line[1] == gene_id:
                modules.append(line)

        if gene_id in expected_modules_constitutive:
            print("Veryify %s constitutive..." % gene_id)
            if expected_modules_constitutive[gene_id]:
                try:
                    assert len(modules) == len(expected_modules_constitutive[gene_id])
                except:
                    print("expt: %d found: %d (%s)" % (len(expected_modules_constitutive[gene_id]), len(modules), gene_id))
                    raise


                for i, mod in enumerate(expected_modules_constitutive[gene_id]):
                    print(modules[i])
                    print(mod)

                    for k, v in mod.items():
                        try:
                            if k == 'lsv_ids':
                                assert all(x in v.split(';') for x in modules[i][expected_headers_constitutive.index(k)].split(';'))
                            else:

                                assert v == modules[i][expected_headers_constitutive.index(k)]
                        except:
                            print("expt: %s found: %s (%d, %s, %s)" % (v, modules[i][expected_headers_constitutive.index(k)], i+1,
                                                                       headers[expected_headers_constitutive.index(k)], gene_id))
                            raise



def verify_mpe(gene_id):
    with open(os.path.join(out_dir, 'mpe_primerable_regions.tsv'), 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')

        headers = next(reader, None)
        mpe_rows = []
        for line in reader:
            if line[1] == gene_id:
                mpe_rows.append(line)

        if gene_id in expected_mpes:
            print("Veryify %s mpe..." % gene_id)
            if expected_mpes[gene_id]:
                try:
                    assert len(mpe_rows) == len(expected_mpes[gene_id])
                except:
                    print("expt: %d found: %d (%s)" % (len(expected_mpes[gene_id]), len(mpe_rows), gene_id))
                    raise


                for i, exectedmperow in enumerate(expected_mpes[gene_id]):
                    print(mpe_rows[i])
                    print(exectedmperow)

                    for expected_header, v in exectedmperow.items():
                        try:
                            assert v == mpe_rows[i][expected_headers_mpe.index(expected_header)]
                        except:
                            print("expt: %s found: %s (%d, %s, %s)" % (v, mpe_rows[i][expected_headers_mpe.index(expected_header)], i+1,
                                                                       headers[expected_headers_mpe.index(expected_header)], gene_id))
                            raise



import sys

def run_tests():

    if len(sys.argv) > 1:
        if sys.argv[-1] in expected_modules:
            run_voila_classify(sys.argv[-1])
            verify_tsvs(sys.argv[-1])
        elif sys.argv[-1] in expected_modules_constitutive:
            run_voila_classify(sys.argv[-1], ['--keep-constitutive'])
            verify_constitutive(sys.argv[-1])

    else:

        run_voila_classify([gene_id for gene_id in expected_modules],
                           additional_args=['--debug'])
        for gene_id in expected_modules:
            verify_tsvs(gene_id)


        run_voila_classify([gene_id for gene_id in expected_modules_constitutive],
                           enabled_outputs="summary,mpe",
                           additional_args=['--keep-constitutive', '--debug'])
        for gene_id in expected_modules_constitutive:
            verify_constitutive(gene_id)
            verify_mpe(gene_id)


    print("Success!")

if __name__ == "__main__":
    run_tests()