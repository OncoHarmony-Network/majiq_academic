import csv
import subprocess
import os, sys

import argparse

parser = argparse.ArgumentParser(description='Majiq long reads test runner')
parser.add_argument('--gene-id', type=str, required=False,
                    help='run test for this specific gene_id instead of all of them')
args = parser.parse_args()

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

def read_result_file(path):
    with open(path, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for line in reader:
            if args.gene_id:
                if args.gene_id == line['gene_id']:
                    yield line
            else:
                yield line


generate_str = "python generate_splice_graph.py -o testcases/sg_generated --json testcases/splice_graphs.json --template-sql ../voila/rna_voila/api/model.sql"
if args.gene_id:
    generate_str += f' -n {args.gene_id}'

subprocess.check_call(generate_str, shell=True)
subprocess.check_call("python main.py --majiq-splicegraph-path testcases/sg_generated.sql --flair-gtf-path testcases/ex.isoforms.gtf --output-path testcases/result -j 1 --debug", shell=True)

errors = False
found_gene = False

for expected, actual in zip(read_result_file('testcases/comparison.tsv'), read_result_file('testcases/result/comparison.tsv')):

    found_gene = True
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
        break

if not found_gene:
    print("Could not find any matching genes!")
elif not errors:
    print("All is well!")
