import cPickle
from scipy import interpolate
import numpy as np
import os
import glob
import ConfigParser

def ConfigSectionMap(Config, section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1



def global_conf_ini(filename):

    global num_experiments, exp_list, readLen, gc_factor, weigh_factor, gc_bins, gc_bins_val, tissue_repl, sam_dir

    Config = ConfigParser.ConfigParser()
    Config.read(filename)
    exp = ConfigSectionMap(Config,"experiments")
    general = ConfigSectionMap(Config,"info")
    exp_list = []
    tissue_repl = {}
    count = 0
    
    print "KKK", exp
    for exp_idx, lstnames in exp.items():
        tissue_repl[exp_idx]=[]
        elist = lstnames.split(',')
        for exp in elist:
            exp_list.append(exp)
            tissue_repl[exp_idx].append(count)
    
    num_experiments = len(exp_list)
    readLen = int(general['readlen'])
    sam_dir = general['samdir']
    gc_factor = [None]*num_experiments
    gc_bins_val = [None]*num_experiments

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
