from majiq.src.config import Config
import majiq.src.filter as majiq_filter
from majiq.grimoire.lsv import SSOURCE, STARGET, InvalidLSV, new_lsv_definition
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage
import majiq.src.utils as majiq_utils
cimport numpy as np
import h5py



cdef bint reliable_in_data(object junc, int exp_idx) except -1:

    cdef object majiq_config = Config()
    cdef int min_read_x_exp = majiq_config.minreads
    cdef int min_npos_x_exp = majiq_config.minpos
    cdef np.ndarray cover = junc.coverage.toarray()[exp_idx]

    return junc.get_read_num(exp_idx) >= min_read_x_exp and np.count_nonzero(cover) >= min_npos_x_exp



cdef int _detect_lsv2(object exon, object gn, str lsv_type, dict dummy, bint only_annot=False) except -1:

    cdef object majiq_config = Config()
    cdef dict sstype = {SSOURCE: ('5prime', 0), STARGET: ('3prime', 1)}
    cdef list jlist
    cdef object x
    cdef object jj, lsv_in
    cdef int counter
    cdef double group_thresh

    jlist = [x for x in exon.get_junctions(sstype[lsv_type][0]) if x is not None and
                                                                   x.get_donor() is not None and
                                                                   x.get_acceptor() is not None]

    if len(jlist) < 2:
        return 0

    for name, ind_list in majiq_config.tissue_repl.items():
        group_thresh = majiq_config.min_exp
        if group_thresh == -1:
            group_thresh = min((len(ind_list) * 0.5), 2)
        counter = 0
        for jj in jlist:
            for exp_idx in ind_list:
                if only_annot or majiq_filter.reliable_in_data(jj, exp_idx,
                                                               minnonzero=majiq_config.minpos,
                                                               min_reads=majiq_config.minreads):
                    counter += 1
            if counter < group_thresh:
                continue
            break
        else:
            continue
        lsv_in = gn.new_lsv_definition(exon, jlist, lsv_type)
        dummy[name][sstype[lsv_type][1]].append(lsv_in)


cdef int _detect_lsv(object exon, str name, str lsv_type, list dummy, bint only_annot=False) except -1:

    cdef object majiq_config = Config()
    cdef dict sstype = {SSOURCE: ('5prime', 0), STARGET: ('3prime', 1)}
    cdef list jlist
    cdef object x
    cdef object jj, lsv_in
    cdef int counter
    cdef double group_thresh

    jlist = [x for x in exon.get_junctions(sstype[lsv_type][0]) if x is not None and
                                                                   x.get_donor() is not None and
                                                                   x.get_acceptor() is not None]

    if len(jlist) < 2:
        return 0
    ind_list = majiq_config.tissue_repl[name]
    group_thresh = majiq_config.min_exp
    if group_thresh == -1:
        group_thresh = min((len(ind_list) * 0.5), 2)

    counter = 0
    for jj in jlist:
        counter = 0
        for exp_idx in ind_list:
            if only_annot or majiq_filter.reliable_in_data(jj, exp_idx,
                                                           minnonzero=majiq_config.minpos,
                                                           min_reads=majiq_config.minreads):
                counter += 1
        if counter < group_thresh:
            continue
        break
    else:
        return 0
    lsv_in = new_lsv_definition(exon, jlist, lsv_type)
    dummy[sstype[lsv_type][1]].append(lsv_in)
    return 1

cdef int __prepare_lsvs(object lsv, str name, object gc_vfunc, list fitfunc_r, object queue) except -1:
    cdef object majiq_config
    cdef int exp_idx
    cdef object gc_f
    cdef dict lsv_q
    cdef object qm

    majiq_config = Config()
    for exp_idx in majiq_config.tissue_repl[name]:
        fname = get_builder_temp_majiq_filename(majiq_config.outDir, majiq_config.sam_list[exp_idx])
        with h5py.File(fname, 'r') as rfa:
            gc_f = None
            if majiq_config.gcnorm:
                try:
                    gc_f = gc_vfunc[exp_idx]
                except:
                    pass

            #majiq_utils.monitor("PRE QUEUE CREAT")
            lsv_q = lsv.to_queue(gc_vfunc=gc_f, fitfunc_r=fitfunc_r[exp_idx], exp=rfa, exp_idx=exp_idx )
            #majiq_utils.monitor("POST QUEUE CREAT")
            qm = QueueMessage(QUEUE_MESSAGE_BUILD_LSV, (lsv_q, exp_idx), 0)
            #qm = QueueMessage(100, (None), 0)
            #queue.put(qm, block=True)

    del lsv



cdef wrap_result_file(object lsv, str name, object gc_vfunc, list fitfunc_r, list lock_per_file):
    cdef object majiq_config
    cdef int exp_idx
    cdef object gc_f

    majiq_config = Config()
    for exp_idx in majiq_config.tissue_repl[name]:
        fname = get_builder_temp_majiq_filename(majiq_config.outDir, majiq_config.sam_list[exp_idx])
        lock_per_file[exp_idx].acquire()
        with h5py.File('%s/%s.boots.hdf5' % (majiq_config.outDir, majiq_config.sam_list[exp_idx]), 'r+') as f:
            gc_f = None
            if majiq_config.gcnorm:
                try:
                    gc_f = gc_vfunc[exp_idx]
                except:
                    pass

            with h5py.File(fname, 'r') as rfa:
                lsv.to_hdf5(hdf5grp=f, gc_vfunc=gc_f, fitfunc_r=fitfunc_r[exp_idx], exp=rfa, exp_idx=exp_idx)

        lock_per_file[exp_idx].release()
    del lsv


cpdef int lsv_detection2(object gn, object gc_vfunc, list fitfunc_r, list lock_per_file, only_real_data=False):
    cdef object majiq_config
    cdef list dummy
    cdef int count
    cdef str name
    cdef list ind_list
    cdef object ex
    cdef object ss, st

    majiq_config = Config()
    # dummy = {}
    count = 0
    for name, ind_list in majiq_config.tissue_repl.items():
        dummy= [[], []]

        for ex in gn.get_exon_list():
            try:
                _detect_lsv(ex, name, SSOURCE, dummy, only_annot=only_real_data)
            except InvalidLSV:
                pass
            try:
                _detect_lsv(ex, name, STARGET, dummy, only_annot=only_real_data)
            except InvalidLSV:
                pass

        for ss in dummy[0]:
            for st in dummy[1]:
                if ss.contained(st):
                    break
            else:
                count += 1
                wrap_result_file(ss, name, gc_vfunc, fitfunc_r, lock_per_file)

        for st in dummy[1]:
            for ss in dummy[0]:
                if st.contained(ss):
                    break
            else:
                count += 1
                wrap_result_file(st, name, gc_vfunc, fitfunc_r, lock_per_file)

    return count


cpdef int lsv_detection(object gn, gc_vfunc, fitfunc_r, queue, only_real_data=False):
    cdef object majiq_config
    cdef list dummy
    cdef int count
    cdef str name
    cdef list ind_list
    cdef object ex
    cdef object ss, st

    majiq_config = Config()
    count = 0
    for name, ind_list in majiq_config.tissue_repl.items():

        dummy = [[], []]

        for ex in gn.get_exon_list():
            try:
                _detect_lsv(ex, name, SSOURCE, dummy, only_annot=only_real_data)
            except InvalidLSV:
                pass
            try:
                _detect_lsv(ex, name, STARGET, dummy, only_annot=only_real_data)
            except InvalidLSV:
                pass

        for ss in dummy[0]:
            for st in dummy[1]:
                if ss.contained(st):
                    break
            else:
                count += 1
                __prepare_lsvs(ss, name, gc_vfunc=gc_vfunc, fitfunc_r=fitfunc_r, queue=queue)

        for st in dummy[1]:
            for ss in dummy[0]:
                if st.contained(ss):
                    break
            else:
                count += 1
                __prepare_lsvs(st, name, gc_vfunc=gc_vfunc, fitfunc_r=fitfunc_r, queue=queue)



    return count

