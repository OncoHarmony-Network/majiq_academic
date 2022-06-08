
import os
import sqlite3
import h5py
import json
import argparse

parser = argparse.ArgumentParser(description='Majiq debug gene generator')
parser.add_argument('-n', '--struct-name', type=str, required=False,
                    help='key in the json file of the struct to create (otherwise all)')
parser.add_argument('-o', '--output-file-prefix', type=str, required=True,
                    help='sql and psi files will start with this')
args = parser.parse_args()

output_splicegraph = args.output_file_prefix + '.sql'
output_psi = args.output_file_prefix + '.psi.voila'

with open('splice_graphs.json', 'r') as f:
    all_structs = json.load(f)

if args.struct_name:
    all_structs = {args.struct_name: all_structs[args.struct_name]}

if os.path.exists(output_splicegraph):
    os.remove(output_splicegraph)
if os.path.exists(output_psi):
    os.remove(output_psi)


con = sqlite3.connect(output_splicegraph)
cur = con.cursor()

with open('voila/rna_voila/api/model.sql') as fp:
    cur.executescript(fp.read())

cur.execute('INSERT INTO genome values (1, "mm10")')
cur.execute('INSERT INTO experiment values ("debugexp")')


f = h5py.File(output_psi, "w")
analysis_type = f.create_dataset("metadata/analysis_type", (), dtype=h5py.special_dtype(vlen=str))
analysis_type[()] = 'psi'
file_version = f.create_dataset("metadata/file_version", (), dtype="<i8")
file_version[()] = 6
experiment_names = f.create_dataset("metadata/experiment_names", (2,1), dtype=h5py.special_dtype(vlen=str))
experiment_names[0, 0] = 'debugexp'
group_names = f.create_dataset("metadata/group_names", (1,), dtype=h5py.special_dtype(vlen=str))
group_names[0] = 'debuggrp'

for struct_name in all_structs:
    struct = all_structs[struct_name]
    exons = struct['exons']
    junctions = struct['junctions']
    introns = struct['introns']
    strand = struct['strand']
    gene_id = struct_name



    cur.execute(f'INSERT INTO gene values ("{gene_id}", "{struct_name}", "{strand}", "Banana")')

    for exon in exons:
        cur.execute(f'INSERT INTO exon values ("{gene_id}", {exon[0]}, {exon[1]}, {exon[0]}, {exon[1]}, 1)')

    for junction in junctions:
        cur.execute(f'INSERT INTO junction values ("{gene_id}", {junction[0]}, {junction[1]}, 1, 1, 0, 1, 1)')
        cur.execute(f'INSERT INTO junction_reads values (10, "debugexp", "{gene_id}", {junction[0]}, {junction[1]})')

    for intron in introns:
        cur.execute(f'INSERT INTO intron_retention values ("{gene_id}", {intron[0]+1}, {intron[1]-1}, 1, 1, 0, 1, 1)')
        cur.execute(f'INSERT INTO intron_retention_reads values (10, "debugexp", "{gene_id}", {intron[0]+1}, {intron[1]-1})')

    f.create_group(f'lsvs/{gene_id}')

con.commit()
con.close()


# lsv_id = 'DEBUG:s:0-0'
# bins = f.create_dataset(f"lsvs/DEBUG/{lsv_id}/bins", (3, 40), dtype="<f4")
# junctions = f.create_dataset(f"lsvs/DEBUG/{lsv_id}/junctions", (3, 2), dtype="<u4")
# lsv_type = f.create_dataset(f"lsvs/DEBUG/{lsv_id}/lsv_type", (), dtype=h5py.special_dtype(vlen=str))
# means = f.create_dataset(f"lsvs/DEBUG/{lsv_id}/means", (3,), dtype="<f4")

f.close()







