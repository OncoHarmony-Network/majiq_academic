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
from pandarallel import pandarallel
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

def extract_junction(row):
    chromStart = row.chromStart
    blockSizes = row.blockSizes.split(',')
    blockStarts = row.blockStarts.split(',')

    temp = [int(chromStart) + int(blockSizes[0])]
    for j in range(len(blockSizes) - 1):
        temp.append(temp[-1] + int(blockStarts[j+1]) - int(blockStarts[j]) - int(blockSizes[j]) + 1)
        temp.append(temp[-1] + int(blockSizes[j+1]) - 1)

    temp = temp[:-1]
    temp = [temp[i:i+2] for i in range(0, len(temp), 2)]

    return temp

def s_parse_bed():


    df = pd.read_csv(args.isq_bed_file, sep='\t', engine='python')#, header=None)
    df = df[df.blockCount != 1]

    pandarallel.initialize(progress_bar=True, verbose=0)
    df['junction'] = df[['chromStart', 'blockSizes', 'blockStarts']].parallel_apply(lambda x: extract_junction(x), axis=1)


    df = df.explode('junction')
    df[['start','end']] = pd.DataFrame(df.junction.tolist(), index= df.index)
    df = df.groupby(['start', 'end']).size().reset_index(name='counts')
    df['junc'] = list(zip(df.start, df.end))
    bed_out = dict(zip(df['junc'], df['counts']))
    return bed_out



def s_parse_gtf():

    df_gtf_all = read_gtf(args.isq_gtf_file).to_pandas()

    df3 = df_gtf_all
    df3 = df3[~df3['feature'].isin(['gene'])]
    df3 = df3.sort_values(['transcript_id','start'], ascending = [False, True]).reset_index(drop=True)
    count = df3['transcript_id'].value_counts()
    df3 = df3[~df3['transcript_id'].isin(count[count < 3].index)]
    df3 = df3[~df3['feature'].isin(['transcript'])]
    new_val = np.nan
    df3.insert(loc=5, column= "next_exon", value=new_val)
    df3['next_exon']= df3.start.shift(-1)
    df3 = df3.drop(df3.groupby('transcript_id').tail(1).index, axis=0)
    df3['next_exon']= df3['next_exon'].astype(int)
    gtf_out = {}
    for i in tqdm(zip(df3['end'], df3['next_exon'])):
        gtf_out[i] = 1

    return gtf_out, df_gtf_all

def s_gtf_bed_combine(bed_out, gtf_out):

    bed_gtf_combine = {}

    for k,v in tqdm(bed_out.items()):
        if gtf_out.get(k):
            bed_gtf_combine[k] = v

    return bed_gtf_combine

def s_parse_tsv():
    df_tsv = pd.read_csv(args.isq_tsv_file, sep='\t', engine='python')
    df_tsv = df_tsv[df_tsv.transcript_id.apply(lambda x: len(x)>1)]
    df_tsv_count = df_tsv.groupby('transcript_id', sort=False).agg({'transcript_id': 'size'}).to_dict()
    return df_tsv_count

def s_final_reads_combine(df_gtf_all, bed_gtf_combine, df_tsv_count):


    transcripts_dict = {}
    junctions_dict = {}

    for gene in tqdm(df_gtf_all.gene_id.unique()):
        df3 = df_gtf_all[df_gtf_all.gene_id == gene]
        df3 = df3[~df3['feature'].isin(['gene'])]
        df3 = df3.sort_values(['transcript_id','start'], ascending = [False, True]).reset_index(drop=True)
        count = df3['transcript_id'].value_counts()
        df3 = df3[~df3['transcript_id'].isin(count[count < 3].index)]
        df3 = df3[~df3['feature'].isin(['transcript'])]
        new_val = np.nan
        df3.insert(loc=5, column= "next_exon", value=new_val)
        df3['next_exon']= df3.start.shift(-1)
        df3 = df3.drop(df3.groupby('transcript_id').tail(1).index, axis=0)

        transcript_dict = {transcript: df_tsv_count['transcript_id'].get(transcript) for transcript in count.index}

        transcripts_dict[gene] = transcript_dict

        junction_dict = {pair: bed_gtf_combine.get(pair) for pair in zip(df3['end'],df3['next_exon'])}
        junctions_dict[gene] = junction_dict

    #out = {'transcript_read_number': transcripts_dict, 'junction_read_number': junctions_dict}
    return transcripts_dict, junctions_dict

        #pickle.dump(out, (open('bed_gtf_tsv_output', 'wb')))


print('~~~Parsing Long Read BED~~~')
bed_out = s_parse_bed()
print('~~~Parsing Long Read GTF~~~')
gtf_out, df_gtf_all = s_parse_gtf()
print('~~~Combining GTF+BED information~~~')
bed_gtf_combine = s_gtf_bed_combine(bed_out, gtf_out)
print('~~~Parsing Long Read TSV~~~')
df_tsv_count = s_parse_tsv()
print('~~~Processing Long Read combined read counts~~~')
transcript_raw_reads, junction_raw_reads = s_final_reads_combine(df_gtf_all, bed_gtf_combine, df_tsv_count)

# transcript_raw_reads format: { 'gene_id': { 'transcript_id' : reads }}
# junction_raw_reads format: { 'gene_id': { (junc_start, junc_end) : reads ) }}

#df_gtf_all = read_gtf(args.isq_gtf_file).to_pandas()
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



    for t_i, (transcript_id, transcript) in enumerate(zip(flairreader.gene_transcript_names(gene_id), flairreader.gene(gene_id))):

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

                for ir_s, ir_e in zip(ir_starts, ir_ends):
                    junc = (ir_s+1, ir_e-1,)
                    transcript_intron_retention.append(junc)
                    transcript_intron_retention_reads.append(junction_raw_reads[gene_id].get(junc, 0))
            else:
                transcript_exons.append((lr_exon[0], lr_exon[1],))

        out_t = {
            'id': transcript_id,
            'strand': strand,
            'experiment': transcript_id,  # f'LR_{gene_id}_{t_i}',
            'exons': transcript_exons,
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
