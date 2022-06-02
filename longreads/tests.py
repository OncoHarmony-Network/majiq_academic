import csv
import subprocess
import os, sys

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

def read_result_file(path):
    with open(path, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for line in reader:
            yield line


subprocess.check_call("python generate_splice_graph.py -o testcases/sg_generated --json testcases/splice_graphs.json --template-sql ../voila/rna_voila/api/model.sql", shell=True)
subprocess.check_call("python main.py --majiq-splicegraph-path testcases/sg_generated.sql --flair-gtf-path /home/pjewell/longreads_debug/collapse3.isoforms.gtf --output-path testcases/result -j 1 --debug", shell=True)

errors = False

for expected, actual in zip(read_result_file('testcases/comparison.tsv'), read_result_file('testcases/result/comparison.tsv')):
    diff = []
    for key in expected.keys():
        if key not in actual:
            diff.append(('Header missing', key, ''))
            continue
        if expected[key] != actual[key]:
            diff.append((key, expected[key], actual[key]))
    if diff:
        errors = True
        print('~~~~~Error in gene:', expected['gene_id'] )
        for key, _exp, _act in diff:
            print('     -', key, 'expected: ', _exp, 'found: ', _act)

if not errors:
    print("All is well!")
