import cPickle
from scipy import interpolate
import numpy as np
import os
import glob

def global_init(readL, my_dir=None, paths=None):

    global num_experiments, exp_list, readLen, gc_factor, weigh_factor, gc_bins, gc_bins_val, tissue_repl

    exp_list = []
    if paths: 
        exp_list = paths
    else:
        print
        print "DIIIIR"
        print
        for path in os.listdir(my_dir):
            print path
            if path.endswith("sam"):
                print "SAAAM", path
                exp_list.append((os.path.join(my_dir, path)))

    print "Experiments:", exp_list

    num_experiments = len(exp_list)
    readLen = readL
    gc_factor = [None]*num_experiments
    gc_bins_val = [None]*num_experiments
