
import pickle
from scipy.sparse import lil_matrix
from grimoire.utils.utils import create_if_not_exists, get_logger
import numpy as np
import random



DELTA_RATIO = 0.2 #TODO to parameters

# UTILS for io
def _find_len(grimoire_obj, logger):
    junc_len = 666 #this value should never be 
    for junction in grimoire_obj:
        if hasattr(junction[0], 'coverage'):
            if type(junction[0].coverage) == lil_matrix:
                event = junction[0].coverage.toarray()
                junc_len = event[0].shape[0]
            else:
                junc_len = junction[0].coverage.shape[0]

            if junc_len != 1: 
                if logger: logger.debug("Junction length is %s"%junc_len)
                return junc_len  

def _find_len2(grimoire_obj, logger):
    junc_len = 666
    for junction in grimoire_obj:
        if hasattr(junction, 'coverage'):
            if type(junction.coverage) == lil_matrix:
                event = junction.coverage.toarray()
                junc_len = event[0].shape[0]
            else:
                junc_len = junction.coverage.shape[0]

            if junc_len != 1: 
                if logger: logger.debug("Junction length is %s"%junc_len)
                return junc_len  

def _load_data_const(grimoire_obj, logger=None):
    CONST_MAX = 2000 # we don't really need more than 2000 constitutive exons 
    """
    Overriding Jordis objects. 
    Should be deleted at some point as the majiq.analysis should instead read them and the object should be extended

    Plus, this should change when the gc_content becomes available again.
    """
    ret = []
    #first get junction length, should be almost always the first entry
    my_len = _find_len(grimoire_obj, logger)

    for i, junction in enumerate(grimoire_obj):
        if i == CONST_MAX: # we don't need more constitutive junctions
            return array(ret)
            
        if type(junction[0].coverage) == lil_matrix:
            event = junction[0].coverage.toarray()
            event = event[0]
        else:
            event = junction[0].coverage

        if event.shape[0] == 1:
            event = [0]*my_len


        ret.append(list(event))

    return array(ret)

def _load_data2(grimoire_obj, logger=None, getnames = False):
    """
    Same as above
    """
    ret = []
    #first get junction length, should be almost always the first entry
    my_len = _find_len2(grimoire_obj, logger)
    for i, junction in enumerate(grimoire_obj):  
        if type(junction.coverage) == lil_matrix:
            event = junction.coverage.toarray()
            event = event[0]
        else:
            event = junction.coverage

        if event.shape[0] == 1:
            event = [0]*my_len


        ret.append(list(event))

    return array(ret)



def _load_data_pair2(data, paired_samples, logger=None):
    "Load the information in data inside paired_samples defaultdict"
    my_len = _find_len(data, logger)
    for i, junction in enumerate(data): #iterate junctions in 1 and add to the inc and exc to the dictionary
        if hasattr(junction[0], 'coverage') and hasattr(junction[1], 'coverage'):

            if type(junction[0].coverage[0]) == lil_matrix:  
                my_inc = junction[0].coverage[0].toarray()
                my_inc = my_inc[0]
            else:
                my_inc = junction[0].coverage[0]

            if type(junction[1].coverage[0]) == lil_matrix:  
                my_exc = junction[1].coverage[0].toarray()
                my_exc = my_exc[0]
            else:
                my_exc = junction[1].coverage[0]

            if my_inc.shape[0] == 1:
                my_inc = [0]*my_len

            if my_exc.shape[0] == 1:
                my_exc = [0]*my_len

            my_name = junction[0].name
            if not my_name:
                my_name = junction[1].name

            paired_samples[my_name].extend([list(my_inc), list(my_exc)])    

def load_data(path, logger=None):
    "Load data from the preprocess step. Could change to a DDBB someday"
    data = pickle.load(open(path))
    return _load_data2(data[1][:,0], logger), _load_data2(data[1][:,1], logger), _load_data_const(data[2], logger) #inc, exc, const

def load_data_n(paths, logger=None, const=False):
    "Load and pair n replicas (for the groups). Doesn't get the constitutive data"
    datas = []
    for path in paths:
        if logger: logger.info("Loading %s..."%path)
        datas.append(pickle.load(open(path))) 

    if logger: logger.info("Combining all matrices into 'paired_samples'...")
    paired_samples = defaultdict(list)
    for data in datas:
        _load_data_pair2(data[1], paired_samples, logger)

    pair_length = len(paths) * 2 # the length of an event that has all replicas is double the number of experiments (each experiment has both inc and exc) 

    #print "PAIRED SAMPLES", len(paired_samples)
    event_names = [] 
    ret = []
    for event_name, junctions in paired_samples.items():
        if len(junctions) == pair_length:
            ret.append(junctions)
            event_names.append(event_name)

    return ret, event_names

def load_data_pair(path1, path2, logger=None, tracklist=None):
    """Pairing functionality should be extracted of this function"""
    #const extracting doesnt change
    data1 = pickle.load(open(path1))
    data2 = pickle.load(open(path2))
    const1 = _load_data_const(data1[2], logger)
    const2 = _load_data_const(data2[2], logger)

    paired_samples = defaultdict(list)
    _load_data_pair2(data1[1], paired_samples, logger)
    _load_data_pair2(data2[1], paired_samples, logger)

    #pair the events: Only keep the paired events
    ret = []
    event_names = [] 
    out_file_borrame = open('%s_%s_all_kirsten.txt'%(os.path.basename(path1), os.path.basename(path2)), 'w')
    for event_name, junctions in paired_samples.items():

        if tracklist:
            if event_name in tracklist:
                logger.info("TRACKLIST (%s): inclusion1: %s \n exclusion1: %s \n inclusion2: %s \n exclusion2: %s \n"%(event_name, junctions[0], junctions[1], junctions[2], junctions[3]))

        if len(junctions) == 4:
            out_file_borrame.write("%s: inc1: %s exc1: %s inc2: %s exc2: %s \n"%(event_name, junctions[0], junctions[1], junctions[2], junctions[3]))
            ret.append(junctions)
            event_names.append(event_name)

    ret = array(ret)
    return ret[:,0], ret[:,1], const1, ret[:,2], ret[:,3], const2, event_names

def load_data_lsv(path, logger=None):
    "Load data from the preprocess step. Could change to a DDBB someday"
    data = pickle.load(open(path))
    lsv_cov_list = []
    lsv_gc = []
    lsv_info = []
    const_list = []
    const_info = []
    const_gc = []
    num_pos = data[1][0].junction_list.shape[1]

    for lsv in data[1]:
        lsv_info.append([lsv.coords,lsv.id, lsv.type,0])
#        cov = np.zeros(shape=(len(lsv.junction_list.shape)), dtype=np.dtype('int'))
#        for ii, lsvcov in enumerate(lsv.junction_list.toarray()):
#            cov[ii,:] = lsvcov
        cov = lsv.junction_list.toarray()
        lsv_cov_list.append( cov )
        gc = lsv.gc_factor.toarray()
        lsv_gc.append(gc)

#    print "LSV COV",lsv_cov_list


    clist = random.sample(data[2], min(5000, len(data[2])))
    const_list = np.zeros(shape=(len(clist),num_pos), dtype=np.dtype('int'))
    const_gc   = np.zeros(shape=(len(clist),num_pos), dtype=np.dtype('float'))
    for cidx,const in enumerate(clist):
        const_info.append(const.id)
        const_list[cidx,:] = const.coverage.toarray()
        const_gc [cidx,:]  = const.gc_factor.toarray()
#        const_list.append(const.coverage.toarray())

    return (lsv_cov_list, lsv_info, lsv_gc), (const_list, const_info, const_gc)

