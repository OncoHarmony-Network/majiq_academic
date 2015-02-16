from scipy import interpolate
import numpy as np
import os
import ConfigParser

global gene_tlb
global gc_factor


def ConfigSectionMap(Config, section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


def keep_info(SEevents, a3, a5, both, SE):
    global SEev, A3SS, A5SS, bothSS, totalSE
    for idx in range(20):
        A3SS[idx] += a3[idx]
        A5SS[idx] += a5[idx]

    for idx, ii in enumerate(SEevents):
        SEev[idx] += ii

    bothSS += both
    totalSE += SE


def print_numbers():
    print "A3SS", A3SS
    print "A5SS", A5SS
    print "SE events", SEev
    print "Total SE", totalSE
    print "BOTH", bothSS


def global_conf_ini(filename, params):

    global num_experiments, exp_list, readLen, tissue_repl, sam_dir, num_mapped_reads, genome, \
        genome_path, outDir, temp_oDir, gene_tlb, strand_specific
    global A3SS, A5SS, SEev, bothSS, totalSE
    global MINREADS, MINPOS

    config = ConfigParser.ConfigParser()
    config.read(filename)
    #TODO: check if filename exists
    exp = ConfigSectionMap(config, "experiments")
    general = ConfigSectionMap(config, "info")
    exp_list = []
    tissue_repl = {}
    temp_oDir = []
    count = 0
    
    MINREADS = params.minreads
    MINPOS = params.minpos
    readLen = int(general['readlen'])
    sam_dir = general['samdir']
    genome = general['genome']
    genome_path = general['genome_path']
    if 'type' in general:
        strand_specific = (general['type'] == 'strand-specific')
    else:
        strand_specific = False
    outDir = params.output
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    for exp_idx, lstnames in exp.items():
        tissue_repl[exp_idx] = []
        elist = lstnames.split(',')
        for exp in elist:
            exp_list.append(exp)
            tissue_repl[exp_idx].append(count)
            count += 1

    num_experiments = len(exp_list)
    num_mapped_reads = [0] * num_experiments
    gene_tlb = {}


    A3SS = [0]*20
    A5SS = [0]*20
    bothSS = 0
    SEev = [0]*5
    totalSE = 0


def get_max_denovo_difference():
    return 1000


def set_gc_factors(bins, factor, means):
    global gc_factor, gc_bins_val, gc_bins, gc_means
    gc_factor = [None]*num_experiments
    for idx, exp in enumerate(exp_list):
#        gc_factor[idx] = interpolate.interp1d( bins[idx], factor[idx],bounds_error=False )
        a = np.append(factor[idx], factor[idx][-1])
        gc_factor[idx] = interpolate.interp1d(means[idx], factor[idx], bounds_error=False, fill_value=1)
#        gc_factor[idx] = interpolate.interp1d( bins[idx], a , bounds_error=False) 
        
    gc_bins_val = factor
    gc_bins = bins
    gc_means = means


def global_init(read_l, my_dir=None, paths=None):

    global num_experiments, exp_list, readLen, gc_factor, gc_bins, gc_bins_val, tissue_repl

    exp_list = []
    if paths: 
        exp_list = paths
    else:
        for path in os.listdir(my_dir):
            print path
            if path.endswith("sam"):
                exp_list.append((os.path.join(my_dir, path)))

    print "Experiments:", exp_list

    num_experiments = len(exp_list)
    readLen = read_l
    gc_factor = [None]*num_experiments
    gc_bins_val = [None]*num_experiments
