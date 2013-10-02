import cPickle
from scipy import interpolate
import numpy as np

num_experiments = 0
exp_list=[]
readLen = 0
gc_factor = []


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
            

    
def global_init(readL):

    global num_experiments, exp_list,readLen, gc_factor,weigh_factor, gc_bins,gc_bins_val, tissue_repl
    exp_list = ['Heart1.chr1','Heart3.chr1','Heart5.chr1','Heart6.chr1',
                'Hippocampus1.chr1','Hippocampus2.chr1','Hippocampus4.chr1','Hippocampus5.chr1','Hippocampus6.chr1',
                'Liver1.chr1','Liver2.chr1','Liver3.chr1','Liver4.chr1','Liver5.chr1','Liver6.chr1',
                'Lung2.chr1','Lung3.chr1','Lung5.chr1','Lung6.chr1',
#                'Spleen1.chr1','Spleen2.chr1','Spleen3.chr1','Spleen5.chr1',
                'Spleen1.chr1','Spleen2.chr1','Spleen5.chr1',
                'Thymus1.chr1','Thymus2.chr1','Thymus3.chr1','Thymus4.chr1','Thymus6.chr1']
 #   exp_list = ['Heart1.chr1','Heart3.chr1',
 #               'Hippocampus1.chr1','Hippocampus2.chr1',
 #               'Liver1.chr1','Liver2.chr1',
 #               'Lung2.chr1','Lung3.chr1',
 #               'Spleen1.chr1','Spleen4.chr1',
 #               'Thymus1.chr1','Thymus2.chr1']
    
#    exp_list =  ['Heart1.chr6','Hippocampus1.chr6']
    #exp_list =  ['Heart1.chr1','Hippocampus1.chr1']
#    exp_list =  ['Hippocampus1.chr1','Liver1.chr1']
#    exp_list =  [ 'Liv_CT22.chr1', 'Liv_CT28.chr1','Liv_CT34.chr1','Liv_CT40.chr1','Liv_CT46.chr1','Liv_CT52.chr1','Liv_CT58.chr1','Liv_CT64.chr1',
#                  'Hyp_CT22.chr1', 'Hyp_CT28.chr1','Hyp_CT34.chr1','Hyp_CT40.chr1','Hyp_CT46.chr1','Hyp_CT52.chr1','Hyp_CT58.chr1','Hyp_CT64.chr1' ]
#    tissue_repl = { 'Liv':[0,1,2,3,4,5,6,7],'Hyp':[8,9,10,11,12,13,14,15] }
#    exp_list = [ 'Hippocampus1.chr1','Hippocampus2.chr1','Hippocampus4.chr1','Hippocampus5.chr1','Hippocampus6.chr1',
#                'Liver1.chr1','Liver2.chr1','Liver3.chr1','Liver4.chr1','Liver5.chr1','Liver6.chr1']
#    tissue_repl = {'Hippocampus.chr1':[0,1,2,3,4],'Liver.chr1':[5,6,7,8,9,10]}
    exp_list =  ['Heart1','Heart3']
    tissue_repl = {'Heart':[0,1,]}
#    exp_list =  ['Heart1.chr1']
    
    weigh_factor = [0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.67221794713393235, 
                    0.822110928783363, 1.0, 0.96555466736662077, 0.97034079788351413, 0.95888448040949026, 0.93916171702266071, 0.9586120159245054, 0.98886869384130682, 0.97031850460351132, 0.97919730108521363, 0.98845338964933505, 0.97762016912650862, 0.96000443463677787, 0.97795543996841916, 0.98483897464546255, 0.98430701054559211, 0.97621404631851538, 0.97557091482162561, 0.99783624218670419, 0.99420256804654417, 0.99996092005852, 1.0, 0.98891003681022327, 0.98408910298925234, 0.98588911669260937, 0.98944197552348001, 0.98861266787997559, 0.98334059809099128, 0.98616818121835459, 0.98356568445706294, 1.0, 1.0, 0.99415876734414588, 1.0, 0.99732413178991319, 0.98571657557568526, 0.98637294512249951, 0.98846297241187242, 1.0, 0.98857076368303576, 1.0, 0.98474007029306976, 0.98212050612598556, 0.99227062085183826, 0.98716235724225032, 0.98604617629365343,
                    0.96908030440229109, 0.97105918006649872, 0.97297718733803484, 0.98431591864639367, 0.98227616224387038, 0.9961571944449884, 0.97565056267585271, 0.96725772937340826, 0.95469906291036666,
                    0.94761567083759468, 0.96284719373281014, 1.0, 0, 0, 0, 0, 0, 0, 0, 0]

    num_experiments = len(exp_list)
    readLen = readL
    gc_factor = [None]*num_experiments
    gc_bins_val = [None]*num_experiments
    with open(r"/home/jordi/working/GCcontent/gc_content_factors.dat", "rb") as input_file:
#    with open(r"/data/ucsc/reads/test_1k/hog/gc_content_factors_Hog.dat", "rb") as input_file:
        temp,gc_bins = cPickle.load(input_file)
#    print temp
    for idx, exp in enumerate(exp_list):
        exp = exp.replace(".chr1","").lower()
#        exp = exp.replace(".chr1","")
#        print temp
        if exp not in temp :
            print "error at global init, not found experiment", exp
            exit(1)
        gc_bins_val[idx] = temp[exp][1]
        a = np.append(temp[exp][1],temp[exp][1][-1])
#        print gc_bins[idx], a
        gc_factor[idx] = interpolate.interp1d(gc_bins[idx], a,bounds_error=False )
