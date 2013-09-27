import cPickle
from scipy import interpolate
import numpy as np
import os


#def global_init(config_file):
#
#
#    global num_experiments, exp_list,readLen, gc_factor,weigh_factor
#   
#
#    fp = open(config_file)
#
#    for ii in fp.readline() :
#        if ii.startswith('#') : continue
#        tab = ii.strip().split('=')
#        if tab[0] == 'EXPERIMENTS' :
#            tab[1].replace('\s','')
            

    
def global_init(readL, path):

    global num_experiments, exp_list,readLen, gc_factor,weigh_factor, gc_bins,gc_bins_val, tissue_repl

    exp_list = []
    for replica_dir in os.listdir(path):
        if replica_dir.find(".sam") == -1: continue
        replica_dir.replace('sorted.sam','')
        exp_list.append(replica_dir)

    print exp_list

    num_experiments = len(exp_list)
    readLen = readL
    gc_factor = [None]*num_experiments
    gc_bins_val = [None]*num_experiments
#    with open(r"/home/jordi/working/GCcontent/gc_content_factors.dat", "rb") as input_file:
##    with open(r"/data/ucsc/reads/test_1k/hog/gc_content_factors_Hog.dat", "rb") as input_file:
#        temp,gc_bins = cPickle.load(input_file)
##    print temp
#    for idx, exp in enumerate(exp_list):
#        exp = exp.replace(".chr1","").lower()
##        exp = exp.replace(".chr1","")
##        print temp
#        if exp not in temp :
#            print "error at global init, not found experiment", exp
#            exit(1)
#        gc_bins_val[idx] = temp[exp][1]
#        a = np.append(temp[exp][1],temp[exp][1][-1])
##        print gc_bins[idx], a
#        gc_factor[idx] = interpolate.interp1d(gc_bins[idx], a,bounds_error=False )
