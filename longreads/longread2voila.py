import argparse
import flairParser
import pickle
import sqlite3
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser(description='Converter from various long read non-majiq output softwares to voila.lr file')
parser.add_argument('-i', '--input-file', type=str, required=True,
                    help='please provide an input text file from a long read program, for example, a gtf')
parser.add_argument('-o', '--output-file', type=str, required=True,
                    help='the path to write the resulting voila file to')
parser.add_argument('-sg', '--splice-graph', type=str, required=True,
                    help='the path to the majiq splice graph file which will be used to align to annotated exons')

args = parser.parse_args()

print('~~~Parsing Flair~~~')
flairreader = flairParser.FlairReader(args.input_file)
print('~~~Done Parsing Flair~~~')

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
print("total genes in LR: ", len(all_gene_ids))
for prog_i, gene_id in enumerate(all_gene_ids):

    # if gene_id != 'ENSG00000256966.6':
    #     continue

    #majiq_gene_id = 'gene:' + gene_id.split('.')[0]
    majiq_gene_id = gene_id
    annotated_exons = IntervalTree.from_tuples((x['start'], x['end'],) for x in sr_gene_exons(majiq_gene_id) if x['start'] != -1 and x['end'] != -1 and x['start'] != x['end'])
    if not annotated_exons:
        continue

    strand = get_strand(majiq_gene_id)

    if prog_i % 100 == 0:
        print(prog_i, len(all_gene_ids), gene_id)
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
