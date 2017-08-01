from majiq.src.config import Config
import majiq.src.filter as majiq_filter
from majiq.grimoire.lsv import SSOURCE, STARGET, InvalidLSV
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage


cdef bint reliable_in_data(object junc, int exp_idx):

    cdef object majiq_config = Config()
    cdef int min_read_x_exp = majiq_config.minreads
    cdef int min_npos_x_exp = majiq_config.minpos
    cdef np.ndarray cover = junc.coverage.toarray()[exp_idx]

    return junc.get_read_num(exp_idx) >= min_read_x_exp and np.count_nonzero(cover) >= min_npos_x_exp



cdef int _detect_lsv(object exon, object gn, str lsv_type, dict dummy, bint only_annot=False):

    cdef object majiq_config = Config()
    cdef dict sstype = {SSOURCE: ('5prime', 0), STARGET: ('3prime', 1)}
    cdef list jlist
    cdef object x
    cdef object jj, lsv_in
    cdef int group_thresh, counter

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


cdef int __prepare_lsvs(object lsv, str name, object gc_vfunc, list fitfunc_r, object queue):
    majiq_config = Config()
    for exp_idx in majiq_config.tissue_repl[name]:
        if majiq_config.gcnorm:
            try:
                gc_f = gc_vfunc[exp_idx]
            except:
                gc_f = None

        lsv_q = lsv.to_queue(gc_vfunc=gc_f, fitfunc_r=fitfunc_r[exp_idx], exp_idx=exp_idx)

        qm = QueueMessage(QUEUE_MESSAGE_BUILD_LSV, (lsv_q, exp_idx), 0)
        queue.put(qm, block=True)


cpdef int lsv_detection(object gn, gc_vfunc, fitfunc_r, queue, only_real_data=False):
    majiq_config = Config()
    dummy = {}
    count = 0
    for name, ind_list in majiq_config.tissue_repl.items():
        dummy[name] = [[], []]

    for ex in gn.get_exon_list():
        try:
            _detect_lsv(ex, gn, SSOURCE, dummy, only_annot=only_real_data)
        except InvalidLSV:
            pass
        try:
            _detect_lsv(ex, gn, STARGET, dummy, only_annot=only_real_data)
        except InvalidLSV:
            pass

    for name, ind_list in majiq_config.tissue_repl.items():

        for ss in dummy[name][0]:
            for st in dummy[name][1]:
                if ss.contained(st):
                    break
            else:
                count += 1
                __prepare_lsvs(ss, name, gc_vfunc=gc_vfunc, fitfunc_r=fitfunc_r, queue=queue)

        for st in dummy[name][1]:
            for ss in dummy[name][0]:
                if st.contained(ss):
                    break
            else:
                count += 1
                __prepare_lsvs(st, name, gc_vfunc=gc_vfunc, fitfunc_r=fitfunc_r, queue=queue)

    return count

