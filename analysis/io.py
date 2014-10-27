import pickle
import numpy as np
import random


def load_data_lsv(path, group_name, logger=None):
    """Load data from the preprocess step. Could change to a DDBB someday"""
    data = pickle.load(open(path))
    lsv_cov_list = []
    lsv_gc = []
    lsv_info = []
    const_info = []
    num_pos = data[1][0].junction_list.shape[1]

    meta_info = data[0]
    meta_info['group'] = group_name
    for lsv in data[1]:
        try:
            lsv_info.append([lsv.coords, lsv.id, lsv.type, 0, lsv.visual])
        except AttributeError, e:
            lsv_info.append([lsv.coords, lsv.id, lsv.type, 0])

        cov = lsv.junction_list.toarray()
        lsv_cov_list.append(cov)
        gc = lsv.gc_factor.toarray()
        lsv_gc.append(gc)

    clist = random.sample(data[2], min(5000, len(data[2])))
    const_list = np.zeros(shape=(len(clist), num_pos), dtype=np.dtype('int'))
    const_gc = np.zeros(shape=(len(clist), num_pos), dtype=np.dtype('float'))
    for cidx, const in enumerate(clist):
        const_info.append(const.id)
        const_list[cidx, :] = const.coverage.toarray()
        const_gc[cidx, :] = const.gc_factor.toarray()

    return meta_info, [lsv_cov_list, lsv_info, lsv_gc], [const_list, const_info, const_gc]
