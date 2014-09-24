#!/usr/bin/env python

from matplotlib import use
use('Agg')

from scipy.stats import beta
import numpy as np
import os
import argparse
import analysis.filter as majiq_filter
import analysis.io as majiq_io
import scipy.stats
import pickle
from random import choice

import ipdb

def res_dump( data, id, path ):

    output = open("%s/%s.naive.pkl"%(path, id), 'w')
    pickle.dump(data, output)

def psi_calc(lsv_junc, args):
    psi = np.zeros(shape=(len(lsv_junc), args.nboots), dtype=np.float)
    for lsv_idx, lsv in enumerate(lsv_junc):

        for iternumber in xrange(args.nboots):
            samples = np.zeros(shape=(len(lsv)),dtype=np.float)
            for i, junction in enumerate( lsv ):
                junction_samples = []
#                ipdb.set_trace()
                for numsamples in xrange(args.npos):
                    junction_samples.append(choice(junction))
                samples[i] = np.sum(junction_samples)

            psi[lsv_idx,iternumber] = samples[0] / np.sum(samples)
    return psi

def calcpsi_func(args):

    print "Calcpsi"
    print "Loading %s..."%args.file
    lsv_junc, const = majiq_io.load_data_lsv(args.file, None) 
    print "Loaded."

    ''' 
        fon[0] = False deactivates num_reads >= 10
        fon[1] = False deactivates npos >= 5
    '''
    fon = [True, False]
    lsv_junc = majiq_filter.lsv_quantifiable( lsv_junc, 1, args.minreads, None , fon)

    psi = psi_calc(lsv_junc[0], args)


    for lsv_idx, lsv in enumerate(psi):
        hi, wgt,nb,b = scipy.stats.histogram(lsv, numbins=nbins)
        ret_psi[lsv_idx] =  hi 

    res_dump([lsv_junc[1],res_psi],id='%s.psi'%args.name,path=args.output)

    return 

def deltapsi_func(args):

    print "Delta psi"
    print "Loading %s..."%args.file1,
    lsv_junc1, const = majiq_io.load_data_lsv(args.file1, None) 
    print "Done."
    print "Loading %s..."%args.file2,
    lsv_junc2, const = majiq_io.load_data_lsv(args.file2, None) 
    print "Done."

    ''' 
        fon[0] = False deactivates num_reads >= 10
        fon[1] = False deactivates npos >= 5
    '''
    fon = [True, False]
    lsv_junc1 = majiq_filter.lsv_quantifiable( lsv_junc1, 1, args.minreads, None , fon)
    lsv_junc2 = majiq_filter.lsv_quantifiable( lsv_junc2, 1, args.minreads, None , fon)

    matched_lsv, matched_info = majiq_filter.lsv_intersection( lsv_junc1, lsv_junc2 )

    psi1 = psi_calc(matched_lsv[0], args)
    psi2 = psi_calc(matched_lsv[1], args)

    delta = np.zeros(shape=(len(matched_info)), dtype=np.float )
    for lsv_idx, info in enumerate(matched_info):
        for itern in xrange(args.nboots):
            dpsi[itern] = psi1[lsv_idx,itern] - psi1[lsv_idx, itern]
        delta[lsv_idx],wgt,nb,b = scipy.stats.histogram(dpsi, numbins=nbins)

    res_dump( [matched_info,delta],id='%s.delta'%args.name,path=args.output )


    return 

if __name__== '__main__':

    parser = argparse.ArgumentParser(description="Script to emulate Naive bootstraping for MAJIQ paper comparition")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('-n', '--num-iterations', default=20000, dest='nboots', type=int, help='Number of iteration of psi for the empirical distribution')
    common.add_argument('-p', '--num-positions', default=50, dest='npos', type=int, help='Number of positions per junction')
    common.add_argument('-r', '--min-reads', default=20, dest='minreads', type=int, help='Minimal number of reads per junction')
    common.add_argument('-i', '--name', default=20, dest='name', type=str, help='Id of the execution')

    calcpsi = argparse.ArgumentParser(add_help=False)
    calcpsi.add_argument('file')
    deltapsi = argparse.ArgumentParser(add_help=False)
    deltapsi.add_argument('file1')
    deltapsi.add_argument('file2')

    subparsers = parser.add_subparsers(help='')
    parser_calcpsi = subparsers.add_parser('calcpsi', help="Calculate PSI values for N experiments, given a folder of preprocessed events by 'majiq preprocess' or SAM/BAM files (This last not implemented yet)", parents=[calcpsi, common])
    parser_calcpsi.set_defaults(func=calcpsi_func)
    parser_deltapair = subparsers.add_parser('deltapsi', help='Calculate Delta PSI values given a pair of experiments (1 VS 1 conditions without replicas)', parents=[common, deltapsi])
    parser_deltapair.set_defaults(func=deltapsi_func)
    args = parser.parse_args()
    args.func(args)



