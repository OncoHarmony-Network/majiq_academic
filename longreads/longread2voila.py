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


def s_parse_bed():
    df = pd.read_csv(args.isq_bed_file, sep='\t', engine='python', header=None)
    df.columns = ['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart',
                  'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']

    # remove 1st row if bed file contained header row
    if df.iloc[0].chromStart == 'chromStart':
        df = df[1:]

    # remove comma from end of blockSizes and blockStarts if existing
    if df.iloc[0].blockSizes[-1] == ',':
        df.blockSizes = df.blockSizes.apply(lambda x: x[:-1])
        df.blockStarts = df.blockStarts.apply(lambda x: x[:-1])

    # iterate over bed row and append start / end values
    junctions = []
    for i in tqdm(range(len(df))):
        chromStart = df.iloc[i]['chromStart']
        blockSizes = df.iloc[i]['blockSizes'].split(',')
        blockStarts = df.iloc[i]['blockStarts'].split(',')

        if len(blockStarts) < 2:
            continue

        temp = [int(chromStart) + int(blockSizes[0])]

        for j in range(len(blockSizes) - 1):
            temp.append(temp[-1] + int(blockStarts[j+1]) - int(blockStarts[j]) - int(blockSizes[j]) + 1)
            temp.append(temp[-1] + int(blockSizes[j+1]) - 1)

        for j in range(len(temp)-1):
            if j % 2 == 0:
                junctions.append(temp[j])
                junctions.append(temp[j+1])

    # count number of start/end pair occurence
    juncs_dict = {}
    for i in range(len(junctions)-1):
        if i % 2 == 0:
            pair = f'{junctions[i]}, {junctions[i+1]}'

            if not juncs_dict.get(pair):
                juncs_dict[pair] = 1
            else:
                juncs_dict[pair] += 1

    # extract all dict key/value/count to a string
    bed_out = '{'
    for key, value in juncs_dict.items():
        bed_out += f'({key}): {value}, '
    bed_out = bed_out[:-2] + '}'

    return bed_out

def s_parse_gtf():

    df_gtf_all = read_gtf(args.isq_gtf_file)

    gtf_out = '{'
    temp = pd.DataFrame(columns=df_gtf_all.columns)

    for transcript_id in tqdm(df_gtf_all.transcript_id.unique()):
        df_gtf = df_gtf_all[df_gtf_all.transcript_id == transcript_id]
        df_gtf.reset_index(inplace=True, drop=True)

        for i in range(len(df_gtf)):

            if df_gtf.iloc[i].feature == 'exon':
                temp.loc[len(temp)] = df_gtf.iloc[i]

            if df_gtf.iloc[i].feature != 'exon':
                if len(temp) > 1:
                    if temp.iloc[0].start < temp.iloc[-1].end:
                        for j in range(len(temp)-1):
                            gtf_out += f'({temp.iloc[j].end}, {temp.iloc[j+1].start}), '
                    else:
                        for j in range(len(temp)-1, 0, -1):
                            gtf_out += f'({temp.iloc[j].end}, {temp.iloc[j-1].start}), '

                temp = pd.DataFrame(columns=df_gtf.columns)

        if len(temp) > 1:
            if temp.iloc[0].start < temp.iloc[-1].end:
                for j in range(len(temp)-1):
                    gtf_out += f'({temp.iloc[j].end}, {temp.iloc[j+1].start}), '
            else:
                for j in range(len(temp)-1, 0, -1):
                    gtf_out += f'({temp.iloc[j].end}, {temp.iloc[j-1].start}), '

        temp = pd.DataFrame(columns=df_gtf.columns)

    gtf_out = gtf_out[:-2] + '}'
    return gtf_out, df_gtf_all

def s_gtf_bed_combine(bed_out, gtf_out):
    bed_out = bed_out[1:-1]
    bed_out = bed_out.replace('(', '')
    bed_out = bed_out.replace('):', ',')
    bed_out = bed_out.replace(' ', '')
    bed_out = bed_out.split(',')

    gtf_out = gtf_out[1:-1]

    bed_in_gtf_out = '{'
    for i in tqdm(range(len(bed_out))):
        if i % 3 == 0:
            pair = f'({bed_out[i]}, {bed_out[i+1]})'

            if pair in gtf_out:
                bed_in_gtf_out += f'{pair}: {bed_out[i+2]}, '


    bed_in_gtf_out = bed_in_gtf_out[:-2] + '}'

    bed_in_gtf_out = bed_in_gtf_out[1:-1]
    bed_in_gtf_out = bed_in_gtf_out.replace('(', '')
    bed_in_gtf_out = bed_in_gtf_out.replace('):', ',')
    bed_in_gtf_out = bed_in_gtf_out.replace(' ', '')
    bed_in_gtf_out = bed_in_gtf_out.split(',')

    final_txt_dict = {}
    for i in tqdm(range(len(bed_in_gtf_out))):
        if i % 3 == 0:
            final_txt_dict[f'({bed_in_gtf_out[i-3]}, {bed_in_gtf_out[i-2]})'] = bed_in_gtf_out[i-1]
    return final_txt_dict

def s_parse_tsv():
    df_tsv = pd.read_csv(args.isq_tsv_file, sep='\t', engine='python')
    df_tsv = df_tsv[df_tsv.transcript_id.apply(lambda x: len(x)>1)]
    return df_tsv

def s_final_reads_combine(df_gtf_all, df_tsv):

    transcript_out = '{\ntranscript_read_number:\n'
    junction_out = 'junction_read_number:\n'

    # for gene in tqdm(df_gtf.gene_id.unique()[0:20]): # outputs initial 20 genes
    for gene in tqdm(df_gtf_all.gene_id.unique()):
        df_gene = df_gtf_all[df_gtf_all.gene_id == gene]

        transcript_out += '{' + gene + ': '
        junction_out += '{' + gene + ': '

        for transcript_id in df_gene.transcript_id.unique():
            if transcript_id == '':
                continue

            transcript_out += f'{transcript_id}: {len(df_tsv[df_tsv.transcript_id == transcript_id])}, '

            df_transcript = df_gene[df_gene.transcript_id == transcript_id]
            df_transcript.reset_index(inplace=True, drop=True)

            temp = pd.DataFrame(columns=df_gene.columns)

            for i in range(len(df_transcript)):

                if df_transcript.iloc[i].feature == 'exon':
                    temp.loc[len(temp)] = df_transcript.iloc[i]

                if df_transcript.iloc[i].feature != 'exon':
                    if len(temp) > 1:
                        if temp.iloc[0].start < temp.iloc[-1].end:
                            for j in range(len(temp)-1):
                                pair = f'({temp.iloc[j].end}, {temp.iloc[j+1].start})'
                                junction_out += f'{pair}: {final_txt_dict.get(pair)}, '
                        else:
                            for j in range(len(temp)-1, 0, -1):
                                pair = f'({temp.iloc[j].end}, {temp.iloc[j-1].start})'
                                junction_out += f'{pair}: {final_txt_dict.get(pair)}, '

                    temp = pd.DataFrame(columns=df_transcript.columns)

            if len(temp) > 1:
                if temp.iloc[0].start < temp.iloc[-1].end:
                    for j in range(len(temp)-1):
                        pair = f'({temp.iloc[j].end}, {temp.iloc[j+1].start})'
                        junction_out += f'{pair}: {final_txt_dict.get(pair)}, '

                else:
                    for j in range(len(temp)-1, 0, -1):
                        pair = f'({temp.iloc[j].end}, {temp.iloc[j-1].start})'
                        junction_out += f'{pair}: {final_txt_dict.get(pair)}, '

            temp = pd.DataFrame(columns=df_transcript.columns)

        transcript_out = transcript_out[:-2] + '},\n'
        junction_out = junction_out[:-2] + '},\n'

    out = transcript_out[:-2] + '\n\n' + junction_out[:-2] + '\n}'
    return out

    #pickle.dump(out, (open('bed_gtf_tsv_output', 'wb')))


print('~~~Parsing Long Read BED~~~')
bed_out = s_parse_bed()
print('~~~Parsing Long Read GTF~~~')
gtf_out, df_gtf_all = s_parse_gtf()
print('~~~Combining GTF+BED information~~~')
final_txt_dict = s_gtf_bed_combine(bed_out, gtf_out)
print('~~~Parsing Long Read TSV~~~')
df_tsv = s_parse_tsv()
print('~~~Processing Long Read final read counts~~~')
reads_out = s_final_reads_combine(df_gtf_all, df_tsv)


flairreader = flairParser.FlairReader(gtf_df=df_gtf_all)


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

def get_strand(gene_id):
    query = conn.execute('SELECT id, name, strand, chromosome FROM gene WHERE id=?', (gene_id,))
    fetch = query.fetchone()
    return fetch[2]


all_genes = {}
all_gene_ids = list(set(flairreader.gene_ids))
#print("total genes in LR: ", len(all_gene_ids))
print('~~~Processing final version of long reads transcript~~~')
for gene_id in tqdm(all_gene_ids):


    #majiq_gene_id = 'gene:' + gene_id.split('.')[0]
    majiq_gene_id = gene_id
    annotated_exons = IntervalTree.from_tuples((x['start'], x['end'],) for x in sr_gene_exons(majiq_gene_id) if x['start'] != -1 and x['end'] != -1 and x['start'] != x['end'])
    if not annotated_exons:
        continue

    strand = get_strand(majiq_gene_id)

    # if prog_i % 100 == 0:
    #     print(prog_i, len(all_gene_ids), gene_id)
    all_genes[gene_id] = []


    for t_i, transcript in enumerate(flairreader.gene(gene_id)):


        transcript_exons = []
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
                transcript_junctions.append((transcript[i-1][1], lr_exon[0]))
                transcript_junctions_reads.append(7)

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

                for ir_s, ir_e in zip(ir_starts, ir_ends):
                    transcript_intron_retention.append((ir_s+1, ir_e-1,))
                    transcript_intron_retention_reads.append(13)
            else:
                transcript_exons.append((lr_exon[0], lr_exon[1],))

        out_t = {
            'strand': strand,
            'experiment': f'LR_{gene_id}_{t_i}',
            'exons': transcript_exons,
            'junctions': transcript_junctions,
            'junction_reads': transcript_junctions_reads,
            'intron_retention': transcript_intron_retention,
            'intron_retention_reads': transcript_intron_retention_reads           
        }

        if strand == '-':
            for key in ('exons', 'junctions', 'junction_reads', 'intron_retention', 'intron_retention_reads',):
                out_t[key].reverse()
        #print(out_t)

        all_genes[gene_id].append(out_t)
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
