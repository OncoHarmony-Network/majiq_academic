import numpy as np

import majiq.grimoire.lsv as majiq_lsv
import majiq.src.config as majiq_config
import majiq.src.filter as majiq_filter
from majiq.grimoire.lsv import SSOURCE, STARGET, InvalidLSV
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage

def reliable_in_data(junc, exp_idx):
    min_read_x_exp = majiq_config.MINREADS
    min_npos_x_exp = majiq_config.MINPOS
    in_data_filter = False
    cover = junc.coverage.toarray()[exp_idx]
    if junc.get_read_num(exp_idx) >= min_read_x_exp and np.count_nonzero(cover) >= min_npos_x_exp:
        in_data_filter = True
    return in_data_filter


def detect_lsv(exon, gn, lsv_type, dummy, only_annot=False):

    sstype = {SSOURCE: ['5prime', 0], STARGET: ['3prime', 1]}
    jlist = exon.get_junctions(sstype[lsv_type][0])
    jlist = [x for x in jlist if x is not None]

    if len(jlist) < 2:
        return
    lsv_in = gn.new_lsv_definition(exon, jlist, lsv_type)

    for name, ind_list in majiq_config.tissue_repl.items():
        group_thresh = majiq_config.min_exp
        if group_thresh == -1:
            group_thresh = min((len(ind_list) * 0.5), 2)
        counter = 0
        for jj in jlist:
            for exp_idx in ind_list:
                if only_annot or majiq_filter.reliable_in_data(jj, exp_idx,
                                                               minnonzero=majiq_config.MINPOS,
                                                               min_reads=majiq_config.MINREADS):
                    counter += 1
            if counter < group_thresh:
                continue
            break
        else:
            continue
        dummy[name][sstype[lsv_type][1]].append(lsv_in)


def wrap_result_queue(lsv, name, gc_vfunc, out_queue, chnk, lsv_list=None, lsv_idx=None):
    qm = QueueMessage(QUEUE_MESSAGE_BUILD_LSV, [majiq_lsv.Queue_Lsv(lsv, name, gc_vfunc), name], chnk)
    out_queue.put(qm, block=True)


def wrap_result_file(lsv, name, gc_vfunc, lsv_list, lsv_idx, lock_per_file=None):
    for dx, exp_idx in enumerate(majiq_config.tissue_repl[name]):

        lock_per_file[exp_idx].acquire()
        lsv_idx[exp_idx] = lsv.to_hdf5(hdf5grp=lsv_list[exp_idx], lsv_idx=lsv_idx[exp_idx], gc_vfunc=gc_vfunc,
                                       exp_idx=dx)
        lsv_list[exp_idx].flush()
        lock_per_file[exp_idx].release()


def lsv_detection(gn, gc_vfunc, lsv_list=None, lsv_idx=None, only_real_data=False, locks=None,
                  logging=None):

    dummy = {}
    for name, ind_list in majiq_config.tissue_repl.items():
        dummy[name] = [[], []]

    for ex in gn.get_exon_list():
        try:
            detect_lsv(ex, gn, SSOURCE, dummy, only_annot=only_real_data)
            detect_lsv(ex, gn, STARGET, dummy, only_annot=only_real_data)
        except InvalidLSV:
            pass

    for name, ind_list in majiq_config.tissue_repl.items():

        for ss in dummy[name][0]:
            for st in dummy[name][1]:
                if ss.contained(st):
                    break
            else:
                wrap_result_file(ss, name, gc_vfunc=gc_vfunc, lsv_list=lsv_list, lsv_idx=lsv_idx, lock_per_file=locks)

        for st in dummy[name][1]:
            for ss in dummy[name][0]:
                if st.contained(ss):
                    break
            else:
                wrap_result_file(st, name, gc_vfunc=gc_vfunc, lsv_list=lsv_list, lsv_idx=lsv_idx, lock_per_file=locks)


