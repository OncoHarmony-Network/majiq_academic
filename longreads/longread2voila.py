import argparse
import flairParser
import pickle
import sqlite3
from intervaltree import Interval, IntervalTree
import pandas as pd
import csv
import sys
import numpy as np
from gtfparse import read_gtf
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='Converter from various long read non-majiq output softwares to voila.lr file')
parser.add_argument('--isq-gtf-file', type=str, required=True,
                    help='For isoquant input, provide the long read GTF file path')
parser.add_argument('--isq-bed-file', type=str, required=True,
                    help='For isoquant input, provide the long read BED file path')
parser.add_argument('--isq-tsv-file', type=str, required=True,
                    help='For isoquant input, provide the long read TSV file path')
parser.add_argument('-o', '--output-file', type=str, required=True,
                    help='the path to write the resulting voila file to')
parser.add_argument('-sg', '--splice-graph', type=str, required=True,
                    help='the path to the majiq splice graph file which will be used to align to annotated exons')
parser.add_argument('--gene-id', type=str, required=False,
                    help='Limit to a gene-id for testing')



args = parser.parse_args()

"""
from collections import namedtuple
_args = namedtuple('args', 'isq_gtf_file isq_bed_file isq_tsv_file output_file splice_graph')
args = _args(
    '/slowdata/longread/raw_isq/00_ENCFF563QZR_sub.transcript_models.gtf',
    '/slowdata/longread/raw_isq/00_ENCFF563QZR_sub.corrected_reads.bed',
    '/slowdata/longread/raw_isq/00_ENCFF563QZR_sub.transcript_model_reads.tsv',
    '/slowdata/longread/ex.isoforms.voila.lr',
    '/slowdata/longread/splicegraph.sql '
)
"""


conn = sqlite3.connect(args.splice_graph)
conn.execute('pragma foreign_keys=ON')

def sr_gene_exons(gene_id):
    query = conn.execute('''
                        SELECT gene_id, start, end, annotated_start, annotated_end, annotated 
                        FROM exon 
                        WHERE gene_id=?
                        ''', (gene_id,))
    while True:
        fetch = query.fetchmany(100)
        if not fetch:
            break
        for x in fetch:
            yield dict(zip(('gene_id', 'start', 'end', 'annotated_start', 'annotated_end', 'annotated'), x))


def reads_new_version(df_gtf, tsv_dict):

    transcripts = {}
    junctions = {}
    exons = {}

    for gene in tqdm(df_gtf['gene_id'].unique()):
        if args.gene_id and gene != args.gene_id:
            continue

        df_gene = df_gtf[df_gtf['gene_id']==gene]
        introns_dict = {}
        junction_read_dict = {}
        exons_read_dict = {}

        junc_pairs_from_sql = [(x['start'], x['end'],) for x in sr_gene_exons(gene) if x['start'] != -1 and x['end'] != -1 and x['start'] != x['end']]
        junc_pairs_from_sql = [(junc_pairs_from_sql[i][1], junc_pairs_from_sql[i+1][0]) for i in range(len(junc_pairs_from_sql)-1)]

        transcripts_list = list(df_gene['transcript_id'].unique())
        transcripts_list.remove('')
        transcripts[gene] = {transcript:tsv_dict.get(transcript) for transcript in sorted(transcripts_list)}

        for transcript in transcripts_list:
            df_transcript = df_gene[df_gene['transcript_id'] == transcript]
            exon_pairs_list = [(row['start'], row['end']) for i,row in df_transcript[1:].iterrows()]

            for junc_pair in junc_pairs_from_sql:
                junc_pair = (junc_pair[0]+1, junc_pair[1]-1)

                for exon_pair in exon_pairs_list:
                    if exon_pair[0] <= junc_pair[0] and junc_pair[1] <= exon_pair[1]:
                        if not introns_dict.get(junc_pair):
                            introns_dict[junc_pair] = tsv_dict.get(transcript)
                        else:
                            introns_dict[junc_pair] += tsv_dict.get(transcript)

            if len(df_transcript) < 3:
                continue

            for i,row in df_transcript.iterrows():
                exon_pair = (row['start'], row['end'])
                if not exons_read_dict.get(exon_pair):
                    exons_read_dict[exon_pair] = tsv_dict.get(transcript)
                else:
                    exons_read_dict[exon_pair] += tsv_dict.get(transcript)

            if df_transcript['strand'].iloc[0] == '-':
                df_transcript['next_exon']= df_transcript.start.shift(1)
                df_transcript = df_transcript[2:]
            else:
                df_transcript['next_exon']= df_transcript.start.shift(-1)
                df_transcript = df_transcript[1:-1]

            for i,row in df_transcript.iterrows():
                pair = (row['end'], int(row['next_exon']))
                if not junction_read_dict.get(pair):
                    junction_read_dict[pair] = tsv_dict.get(transcript)
                else:
                    junction_read_dict[pair] += tsv_dict.get(transcript)


        junctions[gene] = dict(sorted(junction_read_dict.items()))
        junctions[gene].update(dict(sorted(introns_dict.items())))
        exons[gene] = dict(sorted(exons_read_dict.items()))


    return transcripts, junctions, exons

print('~~~Parsing Long Read GTF~~~')
df_gtf = read_gtf(args.isq_gtf_file).to_pandas()

print('~~~Parsing Long Read TSV~~~')
df_tsv = pd.read_csv(args.isq_tsv_file, sep='\t', engine='python')
tsv_dict = df_tsv['transcript_id'].value_counts().to_dict()
print('~~~Processing Long Read combined read counts~~~')
transcript_raw_reads, junction_raw_reads, exons_raw_reads = reads_new_version(df_gtf, tsv_dict)

# transcript_raw_reads format: { 'gene_id': { 'transcript_id' : reads }}
# junction_raw_reads format: { 'gene_id': { (junc_start, junc_end) : reads ) }}

#df_gtf_all = read_gtf(args.isq_gtf_file).to_pandas()
flairreader = flairParser.FlairReader(gtf_df=df_gtf)




def get_strand(gene_id):
    query = conn.execute('SELECT id, name, strand, chromosome FROM gene WHERE id=?', (gene_id,))
    fetch = query.fetchone()
    return fetch[2]


all_genes = {}
all_gene_ids = list(set(flairreader.gene_ids))
#print("total genes in LR: ", len(all_gene_ids))
print('~~~Processing final version of long reads transcript~~~')
for gene_id in tqdm(all_gene_ids):
    if args.gene_id and gene_id != args.gene_id:
        continue


    #majiq_gene_id = 'gene:' + gene_id.split('.')[0]
    majiq_gene_id = gene_id
    annotated_exons = IntervalTree.from_tuples((x['start'], x['end'],) for x in sr_gene_exons(majiq_gene_id) if x['start'] != -1 and x['end'] != -1 and x['start'] != x['end'])
    if not annotated_exons:
        continue

    strand = get_strand(majiq_gene_id)

    # if prog_i % 100 == 0:
    #     print(prog_i, len(all_gene_ids), gene_id)
    all_genes[gene_id] = {'transcripts': [], 'lsvs': {}}



    for t_i, (transcript_id, transcript) in enumerate(zip(flairreader.gene_transcript_names(gene_id), flairreader.gene(gene_id))):

        transcript_exons = []
        transcript_exon_reads = []
        transcript_junctions = []
        transcript_junctions_reads = []
        transcript_intron_retention = []
        transcript_intron_retention_reads = []

        if strand == '-':
            #transcript = [(x[1], x[0]) for x in (reversed(transcript))]
            transcript = [x for x in (reversed(transcript))]


        # detect junctions
        for i, lr_exon in enumerate(transcript):
            if i != 0:
                junc = (transcript[i-1][1], lr_exon[0])
                transcript_junctions.append(junc)
                transcript_junctions_reads.append(junction_raw_reads[gene_id].get(junc, 0))

        # look for exons which completely cross boundry (IR)
        for lr_exon in transcript:
            matching_annotated = annotated_exons.overlap(lr_exon[0], lr_exon[1])
            if len(matching_annotated) > 1:
                # need to break long exon into short exons and IR
                ir_starts = []
                ir_ends = []
                for i, output_exon in enumerate(sorted(matching_annotated)):

                    if i == 0:
                        start = lr_exon[0]
                        end = output_exon.end
                        ir_starts.append(output_exon.end)
                    elif i == len(matching_annotated) - 1:
                        start = output_exon.begin
                        end = lr_exon[1]
                        ir_ends.append(output_exon.begin)
                    else:
                        start = output_exon.begin
                        end = output_exon.end
                        ir_starts.append(output_exon.end)
                        ir_ends.append(output_exon.begin)

                    transcript_exons.append((start, end,))
                    transcript_exon_reads.append(exons_raw_reads[gene_id].get(lr_exon, 0))

                for ir_s, ir_e in zip(ir_starts, ir_ends):
                    junc = (ir_s+1, ir_e-1,)
                    transcript_intron_retention.append(junc)
                    transcript_intron_retention_reads.append(junction_raw_reads[gene_id].get(junc, 0))
            else:
                transcript_exons.append((lr_exon[0], lr_exon[1],))
                # print("transcript", exons_raw_reads[gene_id])
                # print("lr_exon", lr_exon)
                transcript_exon_reads.append(exons_raw_reads[gene_id].get(lr_exon, 0))

        out_t = {
            'id': transcript_id,
            'strand': strand,
            'experiment': transcript_id,  # f'LR_{gene_id}_{t_i}',
            'exons': transcript_exons,
            'exon_reads': transcript_exon_reads,
            'junctions': transcript_junctions,
            'junction_reads': transcript_junctions_reads,
            'intron_retention': transcript_intron_retention,
            'intron_retention_reads': transcript_intron_retention_reads,
            'transcript_reads': transcript_raw_reads[gene_id].get(transcript_id, 0)
        }

        if strand == '-':
            for key in ('exons', 'junctions', 'junction_reads', 'intron_retention', 'intron_retention_reads',):
                out_t[key].reverse()
        #print(out_t)

        all_genes[gene_id]['transcripts'].append(out_t)
        # import pprint
        # pprint.pprint(out_t)
        # break


"""
flair IR detection
for majiq_exon_1, majiq_exon_2 in combinations(zip(all_exons_starts, all_exons_ends), 2):
                # print("ann_start ", annotated_exons_starts)
                # print("ann_end ", annotated_exons_ends)
                # print("1 ",majiq_exon_1[1])
                # print("2 ",majiq_exon_2[0])
                junc_start = majiq_exon_1[1]
                junc_end = majiq_exon_2[0]
                for flair_exon in transcript:
                    # print("flair ", flair_exon)
                    if abs(flair_exon.begin) <= junc_start and abs(flair_exon.end) > junc_end:
                        # print(abs(flair_exon.begin), junc_start, abs(flair_exon.end), junc_end)
                        # print("flair start ",flair_exon.begin)
                        # print("flair end ",flair_exon.end)
                        # print("junc start ",junction_.begin)
                        # print("junc end ",junction_.end)
                        novel_intron = True
                        novel = True
"""

with open(args.output_file, 'wb') as f:
    pickle.dump(all_genes, f)

conn.close()
