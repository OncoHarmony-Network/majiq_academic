import pickle
import numpy as np
import random
from voila.io_voila import VoilaInput
from voila.vlsv import VoilaLsv


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


def dump_lsvs_voila(pickle_path, posterior_matrix, lsvs_info, meta_info, psi_list1=None, psi_list2=None):
    """Create VoilaLSVs objects readable by voila."""
    vlsvs=[]
    psi1, psi2 = None, None
    for ii, bins in enumerate(posterior_matrix):
        lsv_graphic = lsvs_info[ii][-1]
        if psi_list1:
            psi1, psi2 = psi_list1[ii], psi_list2[ii]
        vlsvs.append(VoilaLsv(bins, lsv_graphic=lsv_graphic, psi1=psi1, psi2=psi2))

    pickle.dump(VoilaInput(vlsvs, meta_info), open(pickle_path, 'w'))