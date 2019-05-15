from bisect import bisect_left, bisect_right
from itertools import combinations
from pathlib import Path

from voila import constants
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_utils import generate_means

from voila.classifier.as_types import Graph
from voila.classifier.tsv_writer import TsvWriter
from voila.classifier.tests_expected import expected_modules

from subprocess import Popen, PIPE, STDOUT
import os, shutil
import csv

sg_file = '/home/paul/PycharmProjects/majiq/test_cases/classifier/caleb1/splicegraph.sql'
psi_file = '/home/paul/PycharmProjects/majiq/test_cases/classifier/caleb1/ran_treg.psi.voila'
out_dir = '/home/paul/PycharmProjects/majiq/test_cases/classifier/caleb1/testout'

def run_voila_classify(gene_id):
    os.environ['PYTHONPATH'] = '/home/paul/PycharmProjects/majiq'
    cmd = ('python3',
           '/home/paul/PycharmProjects/majiq/voila/run_voila.py',
           'classify', psi_file, sg_file, '-d', out_dir,
           '--gene-id', gene_id,
           '--decomplexify-threshold', '0.0'
           )

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



def verify_tsvs(gene_id):

    with open(os.path.join(out_dir, 'summary.tsv'), 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')

        headers = next(reader, None)
        modules = []
        for line in reader:
            sum_of_events = sum(int(x) if x else 0 for x in line[2:14])
            if sum_of_events > 1:
                assert line[15] == "True"
            modules.append(line)

        print(modules)
        try:
            assert len(modules) == len(expected_modules[gene_id])
        except:
            print("expt: %d found: %d (%s)" % (len(expected_modules[gene_id]), len(modules), gene_id))
            raise

        for i, mod in enumerate(expected_modules[gene_id]):
            print(mod)

            for j, col in enumerate(mod):
                try:
                    if j == 1:
                        assert all(x in col.split(';') for x in modules[i][j].split(';'))
                    else:
                        assert col == modules[i][j]
                except:
                    print("expt: %s found: %s (%d, %s, %s)" % (col, modules[i][j], i+1, headers[j], gene_id))
                    raise


import sys

def run_tests():

    if len(sys.argv) > 1:
        run_voila_classify(sys.argv[-1])
        verify_tsvs(sys.argv[-1])

    else:
        for gene_id in expected_modules:
            run_voila_classify(gene_id)
            verify_tsvs(gene_id)




    print("Success!")

if __name__ == "__main__":
    run_tests()