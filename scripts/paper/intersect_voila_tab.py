import pickle
import argparse
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='Files to intersect')
parser.add_argument('-t', '--threshold', type=float, default=0.2, dest='thrhold')
args = parser.parse_args()


for ff in args.files:
    fp = open(ff)
    farray = np.loadtxt(ff)
    for ii in farray:
        print ii





