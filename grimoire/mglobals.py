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

    global num_experiments, exp_list, readLen, gc_factor, weigh_factor, gc_bins, gc_bins_val, tissue_repl, sam_dir, num_mapped_reads, genome, outDir, temp_oDir

    Config = ConfigParser.ConfigParser()
    Config.read(filename)
    exp = ConfigSectionMap(Config,"experiments")
    general = ConfigSectionMap(Config,"info")
    exp_list = []
    tissue_repl = {}
    temp_oDir = []
    count = 0
    
    readLen = int(general['readlen'])
    sam_dir = general['samdir']
    outDir = general['output']
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    for exp_idx, lstnames in exp.items():
        tissue_repl[exp_idx]=[]
        elist = lstnames.split(',')
        for exp in elist:
            exp_list.append(exp)
            tissue_repl[exp_idx].append(count)
            temp_oDir.append( "%s/%s"%(outDir,exp))
            if not os.path.exists(temp_oDir[count]):
                os.makedirs(temp_oDir[count])
            count +=1
    num_experiments = len(exp_list)
    

    num_mapped_reads = [0] * num_experiments

    genome = 'mm10'
#    #ad hoc data
#    gc_factor = [None]*num_experiments
#    gc_bins_val = [0]*num_experiments
#
#    with open(r"/home/jordi/working/GCcontent/gc_content_factors.dat", "rb") as input_file:
##    with open(r"/data/ucsc/reads/test_1k/hog/gc_content_factors_Hog.dat", "rb") as input_file:
#        temp,gc_bins = cPickle.load(input_file)
#
#    for idx, exp in enumerate(exp_list):
#        exp = exp.replace(".chr1","").lower()
##        exp = exp.replace(".chr1","")
#        if exp not in temp :
#            print "error at global init, not found experiment", exp
#            exit(1)
#        gc_bins_val[idx] = temp[exp][1]
#        a = np.append(temp[exp][1],temp[exp][1][-1])
##        print gc_bins[idx], a
#        gc_factor[idx] = interpolate.interp1d(gc_bins[idx], a,bounds_error=False )
#
#    weigh_factor = [0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.67221794713393235, 
#                    0.822110928783363, 1.0, 0.96555466736662077, 0.97034079788351413, 0.95888448040949026, 0.93916171702266071, 0.9586120159245054, 0.98886869384130682, 0.97031850460351132, 0.97919730108521363, 0.98845338964933505, 0.97762016912650862, 0.96000443463677787, 0.97795543996841916, 0.98483897464546255, 0.98430701054559211, 0.97621404631851538, 0.97557091482162561, 0.99783624218670419, 0.99420256804654417, 0.99996092005852, 1.0, 0.98891003681022327, 0.98408910298925234, 0.98588911669260937, 0.98944197552348001, 0.98861266787997559, 0.98334059809099128, 0.98616818121835459, 0.98356568445706294, 1.0, 1.0, 0.99415876734414588, 1.0, 0.99732413178991319, 0.98571657557568526, 0.98637294512249951, 0.98846297241187242, 1.0, 0.98857076368303576, 1.0, 0.98474007029306976, 0.98212050612598556, 0.99227062085183826, 0.98716235724225032, 0.98604617629365343,
#                    0.96908030440229109, 0.97105918006649872, 0.97297718733803484, 0.98431591864639367, 0.98227616224387038, 0.9961571944449884, 0.97565056267585271, 0.96725772937340826, 0.95469906291036666,
#                    0.94761567083759468, 0.96284719373281014, 1.0, 0, 0, 0, 0, 0, 0, 0, 0]


def set_gc_factors(bins, factor, means):
    global gc_factor, gc_bins_val, gc_bins, gc_means
    gc_factor = [None]*num_experiments
    for idx, exp in enumerate(exp_list):
#        gc_factor[idx] = interpolate.interp1d( bins[idx], factor[idx],bounds_error=False )
        a = np.append(factor[idx],factor[idx][-1])
        gc_factor[idx] = interpolate.interp1d( means[idx], factor[idx], bounds_error=False) 
#        gc_factor[idx] = interpolate.interp1d( bins[idx], a , bounds_error=False) 
        
    gc_bins_val = factor
    gc_bins = bins
    gc_means = means

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
