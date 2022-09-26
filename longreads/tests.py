import csv
import subprocess
import os, sys

import argparse

parser = argparse.ArgumentParser(description='Majiq long reads test runner')
parser.add_argument('--gene-id', type=str, required=False,
                    help='run test for this specific gene_id instead of all of them')
parser.add_argument('--only-full-gene', action='store_true',
                    help='run tests for only the full length gene comparison (by default, all are run)')
parser.add_argument('--only-modules', action='store_true',
                    help='run tests for only the modules mode gene comparison (by default, all are run)')
# parser.add_argument('--fuzziness5', type=int, default=0,
#                     help='Reduce 5 prime fuzziness of long-read sequencing')
# parser.add_argument('--fuzziness3', type=int, default=0,
#                     help='Reduce 3 prime fuzziness of long-read sequencing')
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

def run_and_check(tmpname, splice_graph_json_path, flair_gtf_path, num5: int, num3: int, ground_truth_path, modules=False):
    generate_str = f"python generate_splice_graph.py -o testcases/{tmpname} --json {splice_graph_json_path} --template-sql ../voila/rna_voila/api/model.sql"
    if args.gene_id:
        generate_str += f' -n {args.gene_id}'

    run_str = f"python main.py --majiq-splicegraph-path testcases/{tmpname}.sql --flair-gtf-path {flair_gtf_path} --fuzziness5 {num5} --fuzziness3 {num3} --output-path testcases/{tmpname}/result -j 1 --debug"
    if modules:
        run_str += ' --per-module'

    subprocess.check_call(generate_str, shell=True)
    subprocess.check_call(run_str, shell=True)

    errors = False
    found_gene = False
    num_success = 0

    for expected, actual in zip(read_result_file(ground_truth_path), read_result_file(f'testcases/{tmpname}/result/comparison.tsv')):

        if expected.get('disabled', '').lower() == 'true':
            continue
        del expected['disabled']
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

        num_success += 1
        # else:
        #     print('Matched!', actual)

    if not found_gene:
        print("Could not find any matching genes!")
    elif not errors:
        print(f"All is well! [{num_success}] successful!")

if not args.only_full_gene:
    run_and_check('testcases_module', 'testcases/testcases_module/splice_graphs.json', 'testcases/testcases_module/ex.isoforms.gtf', 0, 5,'testcases/testcases_module/comparison.tsv')

if not args.only_modules:
    run_and_check('testcases_gene', 'testcases/testcases_gene/splice_graphs.json', 'testcases/testcases_gene/ex.isoforms.gtf', 0, 5, 'testcases/testcases_gene/comparison.tsv')