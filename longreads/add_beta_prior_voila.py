import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import warnings
from scipy.stats import beta
import argparse
import pickle
import h5py
from tqdm import tqdm
from numpy import inf

parser = argparse.ArgumentParser(description='Add beta prior information to long-read voila file, using short read LSV categorization to assign junctions')
parser.add_argument('--lr-voila-file', type=str, required=True,
                    help='Path to long read voila file generated by "longread2voila" script')
parser.add_argument('--lsv-definition-file', type=str, required=True,
                    help='This should be a .psi.voila file which we will match LSV definitions to')
parser.add_argument('--gene-id', type=str, required=False,
                    help='Limit to a gene-id for testing')


args = parser.parse_args()

print('Loading LR data')
with open(args.lr_voila_file, 'rb') as f:
    lr_voila = pickle.load(f)

print('Loading SR data')
sr_voila = h5py.File(args.lsv_definition_file, 'r', driver='core', backing_store=False)

INF_VALUE = 100000

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
        bins[bins == inf] = INF_VALUE
        bins[bins == -inf] = -INF_VALUE

        adj_psi.append(adjusted_psi)
        junc_bins.append(bins.tolist())

        #print(adjusted_psi)
        #ax.plot(x, beta.pdf(x, a, b), 'b-', lw=5, alpha=0.6)

        #plt.savefig('testplot.png')

    return adj_psi, junc_bins


def find_lr_junc_reads(lr_transcripts, junc):
    for lr_transcript in lr_voila[gene_id]['transcripts']:
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

print("Processing LSVS...")
for gene_id in tqdm(sr_voila['lsvs'].keys()):
    if args.gene_id and gene_id != args.gene_id:
        continue

    if gene_id in lr_voila:
        lr_voila[gene_id]['lsvs'] = {}
        for lsv_id in sr_voila['lsvs'][gene_id]:

            lr_reads = []
            njunc = len(sr_voila['lsvs'][gene_id][lsv_id]['junctions'])
            for sr_junc in sr_voila['lsvs'][gene_id][lsv_id]['junctions']:
                sr_junc = tuple(sr_junc)
                reads = find_lr_junc_reads(lr_voila[gene_id]['transcripts'], sr_junc)
                if reads is None:
                    reads = 0
                lr_reads.append(reads)
            if len(lr_reads) != njunc:
                # don't bother quantifying this lsv
                continue

            psi, bins = beta_prior(lr_reads)
            lr_voila[gene_id]['lsvs'][lsv_id] = {'psi': psi, 'bins': bins}

print('writing output')
with open(args.lr_voila_file, 'wb') as f:
    pickle.dump(lr_voila, f)

print("Done")