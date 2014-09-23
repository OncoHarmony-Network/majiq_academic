import pickle
import numpy as np
import random

DELTA_RATIO = 0.2 #TODO to parameters


def load_data_lsv(path, logger=None):
    """Load data from the preprocess step. Could change to a DDBB someday"""
    data = pickle.load(open(path))
    lsv_cov_list = []
    lsv_gc = []
    lsv_info = []
    const_list = []
    const_info = []
    const_gc = []
    num_pos = data[1][0].junction_list.shape[1]

    for lsv in data[1]:
        try:
            lsv_info.append([lsv.coords, lsv.id, lsv.type, 0, lsv.visual])
        except AttributeError, e:
            lsv_info.append([lsv.coords, lsv.id, lsv.type, 0])

#        cov = np.zeros(shape=(len(lsv.junction_list.shape)), dtype=np.dtype('int'))
#        for ii, lsvcov in enumerate(lsv.junction_list.toarray()):
#            cov[ii,:] = lsvcov
        cov = lsv.junction_list.toarray()
        lsv_cov_list.append( cov )
        gc = lsv.gc_factor.toarray()
        lsv_gc.append(gc)

#    print "LSV COV",lsv_cov_list


    clist = random.sample(data[2], min(5000, len(data[2])))
    const_list = np.zeros(shape=(len(clist), num_pos), dtype=np.dtype('int'))
    const_gc = np.zeros(shape=(len(clist), num_pos), dtype=np.dtype('float'))
    for cidx, const in enumerate(clist):
        const_info.append(const.id)
        const_list[cidx, :] = const.coverage.toarray()
        const_gc[cidx, :] = const.gc_factor.toarray()
#        const_list.append(const.coverage.toarray())

    return (lsv_cov_list, lsv_info, lsv_gc), (const_list, const_info, const_gc)

