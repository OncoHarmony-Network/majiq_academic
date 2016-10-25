import numpy as np

import majiq.grimoire.lsv as majiq_lsv
import majiq.src.config as majiq_config
import majiq.src.filter as majiq_filter
from majiq.grimoire.lsv import SSOURCE, STARGET, InvalidLSV
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage
import h5py
import os

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
    jlist = [x for x in jlist if x is not None and x.get_donor() is not None and x.get_acceptor() is not None]
    if len(jlist) < 2:
        return

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
        lsv_in = gn.new_lsv_definition(exon, jlist, lsv_type)
        dummy[name][sstype[lsv_type][1]].append(lsv_in)


def wrap_result_file(lsv, name, gc_vfunc, lsv_list, lock_per_file=None):
    for exp_idx in majiq_config.tissue_repl[name]:

        lock_per_file[exp_idx].acquire()
        f = h5py.File(get_builder_majiq_filename(majiq_config.outDir, lsv_list[exp_idx]),
                      'r+', compression='gzip', compression_opts=9)
        gc_f = None
        if majiq_config.gcnorm:
            gc_f = gc_vfunc[exp_idx]

        f.attrs['data_index'] = lsv.to_hdf5(hdf5grp=f, lsv_idx=f.attrs['data_index'], gc_vfunc=gc_f,
                                exp_idx=exp_idx)
        f.close()
        lock_per_file[exp_idx].release()


def lsv_detection(gn, gc_vfunc, lsv_list=None, only_real_data=False, locks=None, logging=None):

    dummy = {}
    for name, ind_list in majiq_config.tissue_repl.items():
        dummy[name] = [[], []]

    for ex in gn.get_exon_list():
        try:
            detect_lsv(ex, gn, SSOURCE, dummy, only_annot=only_real_data)
        except InvalidLSV:
            pass
        try:
            detect_lsv(ex, gn, STARGET, dummy, only_annot=only_real_data)
        except InvalidLSV:
            pass

    for name, ind_list in majiq_config.tissue_repl.items():

        for ss in dummy[name][0]:
            for st in dummy[name][1]:
                if ss.contained(st):
                    break
            else:
                wrap_result_file(ss, name, gc_vfunc=gc_vfunc, lsv_list=lsv_list, lock_per_file=locks)

        for st in dummy[name][1]:
            for ss in dummy[name][0]:
                if st.contained(ss):
                    break
            else:
                wrap_result_file(st, name, gc_vfunc=gc_vfunc, lsv_list=lsv_list, lock_per_file=locks)


