from rna_voila.config import LongReadsConfig

from rna_voila.voila_log import voila_log
import pickle
import sqlite3
from intervaltree import Interval, IntervalTree
from collections import namedtuple
import pandas as pd
import csv
import sys
import numpy as np
from gtfparse import read_gtf
from tqdm import tqdm
from scipy.stats import beta
from numpy import inf
import warnings
import math
import h5py

"""
Long reads processing
"""

exon = namedtuple('exon', 'start end')
junction = namedtuple('junction', 'start end')
_module = namedtuple('module', 'start end')


def beta_prior(all_junc_reads):

    nj = len(all_junc_reads)
    adj_psi = []
    junc_bins = []

    for junc_i, junc_reads in enumerate(all_junc_reads):
        a = junc_reads + (1/nj)
        b = sum(all_junc_reads[:junc_i] + all_junc_reads[junc_i+1:]) + ((nj-1)/nj) * (nj-1)
        adjusted_psi = a / (a+b)

        #mean, var, skew, kurt = beta.stats(a, b, moments='mvsk')
        #fig, ax = plt.subplots(1, 1)
        x = np.linspace(0,1, 40)
        bins = beta.pdf(x, a, b)

        max_non_inf_value = max(value for value in bins if value != np.inf)
        bins[bins == np.inf] = max_non_inf_value
        bins[bins < 1.0e-20] = 1.0e-20

        first_non_zero = next((x for x in bins if x != 0), None)
        for i in range(len(bins)):
            if bins[i] == 0:
                bins[i] = first_non_zero
            else:
                break

        last_non_zero = next((x for x in reversed(bins) if x != 0), None)
        for i in range(len(bins)-1, -1, -1):
            if bins[i] == 0:
                bins[i] = last_non_zero
            else:
                break

        bins = bins/np.sum(bins)
        adj_psi.append(adjusted_psi)
        junc_bins.append(bins.tolist())

    return adj_psi, junc_bins

def find_lr_junc_reads(lr_transcripts, sr_junc):
    for lr_transcript in lr_transcripts:
        try:
            lr_junc_i = lr_transcript['junctions'].index(sr_junc)
        except ValueError:
            lr_junc_i = None

        try:
            lr_ir_i = lr_transcript['intron_retention'].index(sr_junc)
        except ValueError:
            lr_ir_i = None

        if lr_junc_i is not None and lr_ir_i is not None:
            raise Exception(f"Both IR and Junc found with the same coordinate??, {sr_junc}, {lr_transcript}" )
        elif lr_junc_i is not None:
            return lr_transcript['junction_reads'][lr_junc_i]
        elif lr_ir_i is not None:
            return lr_transcript['intron_retention_reads'][lr_ir_i]

class LRGtfReader:

    def __init__(self, gtf_path=None, gtf_df=None):
        self.modules = {}
        if gtf_path:
            self.df = read_gtf(gtf_path).to_pandas()
        else:
            self.df = gtf_df

    def has_gene(self, gene_id):
        return len(self.df[self.df['gene_id'] == gene_id].index) != 0

    @property
    def gene_ids(self):
        yield from self.df['gene_id']

    def get_exons(self, gene_id, majiq_module_extent=None, modules=False):
        lr_exons = set()
        ord_lr_exons = tuple(x for x in self.gene(gene_id, extent=majiq_module_extent, ignore_starts_ends=True))

        for transcript in ord_lr_exons:
            if modules:
                lr_exons.add(tuple(exon(max(majiq_module_extent[0], e.start) if e.start > 0 else e.start,
                                           min(majiq_module_extent[1], e.end) if e.end > 0 else e.end) for e in
                                      transcript))
            else:
                lr_exons.add(tuple(exon(e.start, e.end) for e in transcript))

        return lr_exons

    def extend_modules(self, module_extents, lr_transcripts):

        done_processing = False
        while not done_processing:

            for transcript in lr_transcripts:
                for i in range(len(transcript) - 1):
                    junc = junction(transcript[i].end, transcript[i + 1].start)

                    start_in_module = None
                    end_in_module = None
                    for module in module_extents:
                        if junc.start >= module.start and junc.start <= module.end:
                            start_in_module = module
                        if junc.end >= module.start and junc.end <= module.end:
                            end_in_module = module

                    if start_in_module and end_in_module and start_in_module != end_in_module:
                        module_extents.remove(start_in_module)
                        module_extents.remove(end_in_module)
                        module_extents.append(_module(start_in_module.start, end_in_module.end))
                        break

                    if start_in_module and not end_in_module:
                        # extend to the next exon
                        module_extents.remove(start_in_module)
                        module_extents.append(_module(start_in_module.start, junc.end))
                        break

                    if not start_in_module and end_in_module:
                        # extend to the previous exon
                        module_extents.remove(end_in_module)
                        module_extents.append(_module(junc.start, end_in_module.end))
                        break

                else:
                    continue
                break
            else:
                done_processing = True

        return module_extents

    def gene_transcript_names(self, gene_id):
        df_gene = self.df[self.df['gene_id'] == gene_id]
        for index, row in df_gene.iterrows():
            if row.feature == 'transcript':
                yield row['transcript_id']

    def gene(self, gene_id, extent=None, ignore_starts_ends=False):

        df_gene = self.df[self.df['gene_id'] == gene_id]

        found_transcripts = set()
        transcript_exons = []

        def append_next(_transcript_exons):
            if _transcript_exons:
                # _transcript_exons = sorted(_transcript_exons, key=lambda e: e.start)
                if ignore_starts_ends:
                    _transcript_exons[0] = exon(-_transcript_exons[0].start, _transcript_exons[0].end)
                    _transcript_exons[-1] = exon(_transcript_exons[-1].start, -_transcript_exons[-1].end)
                _transcript_exons = tuple(_transcript_exons)
                if not _transcript_exons in found_transcripts:
                    found_transcripts.add(_transcript_exons)
                    return tuple(_transcript_exons)

        for index, row in df_gene.iterrows():

            if row.feature == 'transcript':
                ret = append_next(transcript_exons)
                if ret:
                    yield ret

                transcript_exons = []
                continue

            elif row.feature == 'exon':
                _exon = exon(int(row.start), int(row.end))

                if extent:
                    if (_exon.end >= extent[0] and _exon.end <= extent[1]) or \
                            (_exon.start >= extent[0] and _exon.start <= extent[1]):
                        transcript_exons.append(_exon)
                else:
                    transcript_exons.append(_exon)

        ret = append_next(transcript_exons)  # catch the last append
        if ret:
            yield ret

def longReadsInputsToLongReadsVoila():
    config = LongReadsConfig()
    log = voila_log()

    log.info("Running LR")

    if not config.only_update_psi:

        conn = sqlite3.connect(config.splice_graph_file)
        conn.execute('pragma foreign_keys=ON')

        import warnings
        warnings.filterwarnings('ignore')

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
                if config.gene_id and gene != config.gene_id:
                    continue

                df_gene = df_gtf[df_gtf['gene_id'] == gene]
                introns_dict = {}
                junction_read_dict = {}
                exons_read_dict = {}

                junc_pairs_from_sql = [(x['start'], x['end'],) for x in sr_gene_exons(gene) if
                                       x['start'] != -1 and x['end'] != -1 and x['start'] != x['end']]
                junc_pairs_from_sql = [(junc_pairs_from_sql[i][1], junc_pairs_from_sql[i + 1][0]) for i in
                                       range(len(junc_pairs_from_sql) - 1)]

                transcripts_list = list(df_gene['transcript_id'].unique())
                transcripts_list.remove('')
                transcripts[gene] = {transcript: tsv_dict.get(transcript) for transcript in sorted(transcripts_list)}

                for transcript in transcripts_list:
                    df_transcript = df_gene[df_gene['transcript_id'] == transcript]
                    exon_pairs_list = [(row['start'], row['end']) for i, row in df_transcript[1:].iterrows()]

                    for junc_pair in junc_pairs_from_sql:
                        junc_pair = (junc_pair[0] + 1, junc_pair[1] - 1)

                        for exon_pair in exon_pairs_list:
                            if exon_pair[0] <= junc_pair[0] and junc_pair[1] <= exon_pair[1]:
                                if not introns_dict.get(junc_pair):
                                    introns_dict[junc_pair] = tsv_dict.get(transcript)
                                else:
                                    introns_dict[junc_pair] += tsv_dict.get(transcript)

<<<<<<< HEAD
                    if len(df_transcript) < 3:
                        continue
=======
    log.info('~~~Parsing Long Read TSV~~~')
    df_tsv = pd.read_csv(config.lr_tsv_file, sep='\t', engine='python')
    df_tsv.columns.values[0] = 'transcript_id'
    df_tsv.columns.values[1] = 'count'
    df_tsv['count'] = df_tsv['count'].apply(lambda x: math.ceil(x))
    df_tsv.set_index('transcript_id', inplace=True)
    tsv_dict = df_tsv['count'].to_dict()
    log.info('~~~Processing Long Read combined read counts~~~')
    transcript_raw_reads, junction_raw_reads, exons_raw_reads = reads_new_version(df_gtf, tsv_dict)
>>>>>>> 3bd6927d (transcript read file format change)

                    for i, row in df_transcript.iterrows():
                        exon_pair = (row['start'], row['end'])
                        if not exons_read_dict.get(exon_pair):
                            exons_read_dict[exon_pair] = tsv_dict.get(transcript)
                        else:
                            exons_read_dict[exon_pair] += tsv_dict.get(transcript)

                    if df_transcript['strand'].iloc[0] == '-':
                        df_transcript['next_exon'] = df_transcript.start.shift(1)
                        df_transcript = df_transcript[2:]
                    else:
                        df_transcript['next_exon'] = df_transcript.start.shift(-1)
                        df_transcript = df_transcript[1:-1]

                    for i, row in df_transcript.iterrows():
                        pair = (row['end'], int(row['next_exon']))
                        if not junction_read_dict.get(pair):
                            junction_read_dict[pair] = tsv_dict.get(transcript)
                        else:
                            junction_read_dict[pair] += tsv_dict.get(transcript)

                junctions[gene] = dict(sorted(junction_read_dict.items()))
                junctions[gene].update(dict(sorted(introns_dict.items())))
                exons[gene] = dict(sorted(exons_read_dict.items()))

            return transcripts, junctions, exons

        log.info('~~~Parsing Long Read GTF~~~')
        df_gtf = read_gtf(config.lr_gtf_file).to_pandas()

        log.info('~~~Parsing Long Read TSV~~~')
        df_tsv = pd.read_csv(config.lr_tsv_file, sep='\t', engine='python')
        tsv_dict = df_tsv['transcript_id'].value_counts().to_dict()
        log.info('~~~Processing Long Read combined read counts~~~')
        transcript_raw_reads, junction_raw_reads, exons_raw_reads = reads_new_version(df_gtf, tsv_dict)

        # transcript_raw_reads format: { 'gene_id': { 'transcript_id' : reads }}
        # junction_raw_reads format: { 'gene_id': { (junc_start, junc_end) : reads ) }}

        # df_gtf_all = read_gtf(args.isq_gtf_file).to_pandas()
        lrreader = LRGtfReader(gtf_df=df_gtf)

        def get_strand(gene_id):
            query = conn.execute('SELECT id, name, strand, chromosome FROM gene WHERE id=?', (gene_id,))
            fetch = query.fetchone()
            return fetch[2]

        all_genes = {}
        all_gene_ids = list(set(lrreader.gene_ids))

        log.info('~~~Processing final version of long reads transcript~~~')
        for gene_id in tqdm(all_gene_ids):
            if config.gene_id and gene_id != config.gene_id:
                continue

            # majiq_gene_id = 'gene:' + gene_id.split('.')[0]
            majiq_gene_id = gene_id
            annotated_exons = IntervalTree.from_tuples((x['start'], x['end'],) for x in sr_gene_exons(majiq_gene_id) if
                                                       x['start'] != -1 and x['end'] != -1 and x['start'] != x['end'])
            if not annotated_exons:
                continue

            strand = get_strand(majiq_gene_id)

            all_genes[gene_id] = {'transcripts': [], 'lsvs': {}}

            for t_i, (transcript_id, transcript) in enumerate(
                    zip(lrreader.gene_transcript_names(gene_id), lrreader.gene(gene_id))):

                transcript_exons = []
                transcript_exon_reads = []
                transcript_junctions = []
                transcript_junctions_reads = []
                transcript_intron_retention = []
                transcript_intron_retention_reads = []

                if strand == '-':
                    # transcript = [(x[1], x[0]) for x in (reversed(transcript))]
                    transcript = [x for x in (reversed(transcript))]

                # detect junctions
                for i, lr_exon in enumerate(transcript):
                    if i != 0:
                        junc = (transcript[i - 1][1], lr_exon[0])
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
                            junc = (ir_s + 1, ir_e - 1,)
                            transcript_intron_retention.append(junc)
                            transcript_intron_retention_reads.append(junction_raw_reads[gene_id].get(junc, 0))
                    else:
                        transcript_exons.append((lr_exon[0], lr_exon[1],))
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

                all_genes[gene_id]['transcripts'].append(out_t)

        conn.close()

    else:
        log.info("Deleting previous Beta Priors...")
        with open(config.output_file, 'rb') as f:
            all_genes = pickle.load(f)
            for gene_id in all_genes:
                del all_genes[gene_id]['lsvs']

    if config.voila_file:
        sr_voila = h5py.File(config.voila_file, 'r', driver='core', backing_store=False)
        log.info("Processing Beta Priors...")
        for gene_id in tqdm(sr_voila['lsvs'].keys()):
            if config.gene_id and gene_id != config.gene_id:
                continue

            if gene_id in all_genes:
                all_genes[gene_id]['lsvs'] = {}
                for lsv_id in sr_voila['lsvs'][gene_id]:

                    lr_reads = []

                    for sr_junc in sr_voila['lsvs'][gene_id][lsv_id]['junctions']:
                        sr_junc = tuple(sr_junc)
                        reads = find_lr_junc_reads(all_genes[gene_id]['transcripts'], sr_junc)
                        if reads is None:
                            reads = 0
                        lr_reads.append(reads)

                    psi, bins = beta_prior(lr_reads)
                    all_genes[gene_id]['lsvs'][lsv_id] = {'psi': psi, 'bins': bins}

    with open(config.output_file, 'wb') as f:
        pickle.dump(all_genes, f)



