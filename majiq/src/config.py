import os
import ConfigParser
from scipy import interpolate
import numpy as np


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


def global_conf_ini(filename, params, only_db=False):
    global num_experiments, exp_list, readLen, tissue_repl, sam_dir, num_mapped_reads, genome, \
        genome_path, outDir, temp_oDir, gene_tlb, strand_specific, permissive_ir, gcnorm, dbfile
    global A3SS, A5SS, SEev, bothSS, totalSE
    global MINREADS, MINPOS, MIN_INTRON
    global num_final_chunks, min_denovo, nrandom_junctions

    if not only_db:
        num_final_chunks = params.nthreads  if params.nthreads > 1 else 1
    else:
        num_final_chunks = 1
    min_denovo = params.min_denovo
    gcnorm = params.gcnorm
    nrandom_junctions = 5000
    config = ConfigParser.ConfigParser()
    config.read(filename)
    # TODO: check if filename exists
    exp = ConfigSectionMap(config, "experiments")
#    lengths_exp = ConfigSectionMap(config, "readlen")
    general = ConfigSectionMap(config, "info")
    exp_list = []
    tissue_repl = {}

    temp_oDir = []
    count = 0

    MINREADS = params.minreads
    MINPOS = params.minpos

    if not only_db:
        permissive_ir = params.permissive
        MIN_INTRON = params.min_intronic_cov

    #readLen = [int(xx) for xx in general['readlen'].split([','])
    sam_dir = general['samdir']
    genome = general['genome']
    genome_path = general['genome_path']
    readLen = int(general['readlen'])


    if 'type' in general:
        strand_specific = (general['type'] == 'strand-specific')
    else:
        strand_specific = False
    outDir = params.output
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    dbfile = "%s/tmp/db.hdf5" % outDir

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

    sam_list = []
    for exp_idx, exp in enumerate(exp_list):
        samfile = "%s/%s.bam" % (sam_dir, exp)
        if not os.path.exists(samfile):
            raise RuntimeError("Skipping %s.... not found" % samfile)
        baifile = "%s/%s.bam.bai" % (sam_dir, exp)
        if not os.path.exists(baifile):
            raise RuntimeError("Skipping %s.... not found ( index file for bam file is required)" % baifile)
        sam_list.append(samfile)
        exp_list[exp_idx] = os.path.split(exp)[1]

    return sam_list


def global_default():
    global num_experiments, exp_list, readLen, tissue_repl, sam_dir, num_mapped_reads, genome, \
        genome_path, outDir, temp_oDir, gene_tlb, strand_specific, permissive_ir, gc_norm
    global A3SS, A5SS, SEev, bothSS, totalSE
    global MINREADS, MINPOS, MIN_INTRON

    num_experiments = 1
    readLen = 30
    MINREADS = 1
    MINPOS = 1
    genome = 'mm10'

    gene_tlb = {}


def get_max_denovo_difference():
    return 500


def add_chunk():
    global num_final_chunks
    num_final_chunks += 1

##TODO: Deprecated
def set_gc_factors(bins, factor, means):
    global gc_factor, gc_bins_val, gc_bins, gc_means
    gc_factor = [None] * num_experiments
    for idx, exp in enumerate(exp_list):
        a = np.append(factor[idx], factor[idx][-1])
        gc_factor[idx] = interpolate.interp1d(means[idx], factor[idx], bounds_error=False, fill_value=1)

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
    gc_factor = [None] * num_experiments
    gc_bins_val = [None] * num_experiments
