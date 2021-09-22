from bisect import bisect_left, bisect_right
from itertools import combinations
from pathlib import Path

from rna_voila import constants
from rna_voila.api import SpliceGraph, Matrix
from rna_voila.api.matrix_utils import generate_means

from rna_voila.classifier.as_types import Graph
from rna_voila.classifier.tsv_writer import TsvWriter


from subprocess import Popen, PIPE, STDOUT
import os, shutil
import csv



# this should vary depending on the group of tests to run
# changed relative import path here:
try:
    from .tests_expected_t_cells_2 import *
except:
    from tests_expected_t_cells_2 import *
import tempfile
#from tests_expected_t_cells_3 import *


if 'VOILA_TEST_OUTPUT_DIR' in os.environ:
    out_dir = os.environ['VOILA_TEST_OUTPUT_DIR']
else:
    out_dir = tempfile.mkdtemp()

os.makedirs(out_dir, exist_ok=True)

majiq_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

def run_voila_tsv(gene_ids, additional_args=[]):
    if os.environ.get('JENKINS_HOME', None):
        cmd = ['voila']
    else:
        os.environ['PYTHONPATH'] = majiq_dir
        cmd = ['python3', os.path.join(majiq_dir, 'rna_voila', 'run_voila.py')]
    cmd += ['tsv', '-f', os.path.join(out_dir, 'comp_tsv.tsv'), psi_file, sg_file]

    for arg in additional_args:
        cmd.append(arg)
    cmd.append('--gene-ids')
    cmd += gene_ids

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

def read_voila_tsv_juncs():
    with open(os.path.join(out_dir, 'comp_tsv.tsv'), 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')
        headers = next(reader, None)
        for line in reader:
            yield line


def read_voila_classify_juncs():
    with open(os.path.join(out_dir, 'junctions.tsv'), 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')
        headers = next(reader, None)
        for line in reader:
            yield line


def run_voila_classify(gene_ids, enabled_outputs='all', additional_args=[]):
    if os.environ.get('JENKINS_HOME', None):
        cmd = ['voila']
    else:
        os.environ['PYTHONPATH'] = majiq_dir
        cmd = ['python3', os.path.join(majiq_dir, 'rna_voila', 'run_voila.py')]
    cmd += ['modulize', psi_file, sg_file, '-d', out_dir,
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


expected_headers = ['module_id', 'gene_id', 'gene_name', 'seqid','strand', 'lsv_ids', 'cassette_exon', 'tandem_cassette', 'alt3ss', 'alt5ss', 'p_alt3ss',
                    'p_alt5ss', 'alt3and5ss', 'mutually_exclusive', 'alternative_intron', 'ale', 'afe', 'p_ale', 'p_afe', 'orphan_junction', 'other',
                    'multi_exon_spanning',   'complex', 'number-of-events']
expected_headers_constitutive = ['module_id', 'gene_id', 'gene_name', 'seqid','strand', 'lsv_ids', 'cassette_exon', 'tandem_cassette', 'alt3ss', 'alt5ss', 'p_alt3ss',
                    'p_alt5ss', 'alt3and5ss', 'mutually_exclusive', 'alternative_intron', 'ale', 'afe', 'p_ale', 'p_afe', 'orphan_junction', 'other',
                                 'constitutive_junction', 'constitutive_intron',
                    'multi_exon_spanning',   'complex', 'number-of-events']
expected_headers_mpe= ['module_id', 'gene_id', 'gene_name', 'seqid','strand', 'lsv_id',
                       'module_event_combination',"Type","Edge of the Module",
                       "Reference Exon Coord","Reference Exon De Novo","Reference Exon Exitrons","Reference Exon Constant Region",
                       "Reference Exon Trimmed","Constitutive Direction","Constitutive Regions",
                       "Constitutive De Novo","Constitutive Exon or Intron"]


# map of globals from the expected test cases file to tsv filenames to read
quant_verif_groups = {'expected_alternative_intron': 'alternative_intron.tsv',
                      'expected_cassette_exons': 'cassette.tsv',
                      'expected_alt3ss': 'alt3prime.tsv',
                      'expected_alt5ss': 'alt5prime.tsv',}


def verify_tsvs(gene_id):

    for quant_verification in quant_verif_groups:
        # check if the group of tests is defined in the expected file
        if quant_verification in globals():
            quants_to_verify = globals().get(quant_verification)


            with open(os.path.join(out_dir, quant_verif_groups[quant_verification]), 'r', newline='') as csvfile:
                reader = csv.reader(csvfile, dialect='excel-tab', delimiter='\t')
                headers = next(reader, None)
                events = []
                for line in reader:
                    if line[1] == gene_id:
                        events.append(line)


                if gene_id in quants_to_verify:
                    for i, mod in enumerate(quants_to_verify[gene_id]):
                        print(mod)

                        for k, v in mod.items():
                            try:
                                assert v == events[i][headers.index(k)]
                            except:
                                print("expt: %s found: %s (%d, %s, %s)" % (v, events[i][headers.index(k)], i+1,
                                                                           events[i][headers.index(k)], gene_id))
                                raise


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

def test_run_functional_tests():
    #print('test run', expected_modules, quant_verif_groups)

    retain_output = '-k' in sys.argv
    if retain_output:
        sys.argv.remove('-k')
    try:
        if len(sys.argv) > 1:
            if sys.argv[-1] in expected_modules:
                run_voila_classify(sys.argv[-1])
                verify_tsvs(sys.argv[-1])
            elif sys.argv[-1] in expected_modules_constitutive:
                run_voila_classify(sys.argv[-1], ['--keep-constitutive'])
                verify_constitutive(sys.argv[-1])

        else:




            gene_ids = list(expected_modules.keys())
            print(gene_ids)
            for quant_verification in quant_verif_groups:
                # check if the group of tests is defined in the expected file
                if quant_verification in globals():
                    gene_ids += list(globals().get(quant_verification).keys())

            # run_voila_tsv(gene_ids,
            #               additional_args=['--debug'])
            # lsv_ids_from_voila_tsv = set()
            # lsv_ids_from_voila_classify = set()
            # for l in read_voila_tsv_juncs():
            #     lsv_ids_from_voila_tsv.add(l[2])
            #
            #
            # run_voila_classify(gene_ids,
            #                    additional_args=['--debug', '--output-complex'])
            # for l in read_voila_classify_juncs():
            #     lsv_ids_from_voila_classify.add(l[5])
            #
            # print(len(lsv_ids_from_voila_classify))
            # print(len(lsv_ids_from_voila_tsv))
            #
            # return

            run_voila_classify(gene_ids,
                               additional_args=['--debug'])
            for gene_id in gene_ids:
                verify_tsvs(gene_id)


            return


            run_voila_classify([gene_id for gene_id in expected_modules_constitutive],
                               enabled_outputs="summary,mpe,events,junctions",
                               additional_args=['--keep-constitutive',
                                                '--debug',
                                                '--output-complex',
                                                '--keep-no-lsv'])
            for gene_id in expected_modules_constitutive:
                verify_constitutive(gene_id)
                verify_mpe(gene_id)

        print("Success!")
        exit_code = 0
    except:
        print("Some test failed!")
        import traceback
        print(traceback.format_exc())
        exit_code = 1
    finally:
        if retain_output:
            print("Output folder %s was retained" % out_dir)
        else:
            shutil.rmtree(out_dir)

    sys.exit(exit_code)

if __name__ == "__main__":
    test_run_functional_tests()