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


def detect_lsv(exon, gn, lsv_type, dummy, jun, only_annot=False):

    sstype = {SSOURCE: ['5prime', 0], STARGET: ['3prime', 1]}
    jlist = exon.get_junctions(sstype[lsv_type][0])
    jlist = [x for x in jlist if x is not None]
    if jlist < 2:
        return
    lsv_in = gn.new_lsv_definition(exon, jlist, lsv_type)

    for name, ind_list in majiq_config.tissue_repl.items():
        group_thresh = majiq_config.min_exp
        if group_thresh == -1:
            group_thresh = min((len(ind_list) * 0.5), 2)
        counter = 0
        e_data = 0
        for jj in jlist:
            for exp_idx in ind_list:
                if only_annot or majiq_filter.reliable_in_data(jj, exp_idx,
                                                               minnonzero=majiq_config.MINPOS,
                                                               min_reads=majiq_config.MINREADS):
                    counter += 1
            if counter < group_thresh:
                continue
            e_data += 1
            try:
                jun[name].add(jj)
            except KeyError:
                jun[name] = set()
                jun[name].add(jj)
        if e_data == 0:
            continue
        dummy[name][sstype[lsv_type][1]].append(lsv_in)
        return


def lsv_detection(gn, gc_vfunc, chnk, only_real_data=False, out_queue=None, logging=None):

    const_set = {}

    local_const = set(gn.get_all_junctions())
    local_lsv_jun = {}

    dummy = {}
    for name, ind_list in majiq_config.tissue_repl.items():
        dummy[name] = [[], []]
        const_set[name] = set()

    for ex in gn.get_exon_list():
        try:
            detect_lsv(ex, gn, SSOURCE, dummy, local_lsv_jun)
        except InvalidLSV:
            pass

        try:
            detect_lsv(ex, gn, STARGET, dummy, local_lsv_jun)
        except InvalidLSV:
            pass

    for name, ind_list in majiq_config.tissue_repl.items():

        local_const.difference(local_lsv_jun)
        const_set[name].update(local_const)

        njuncs = len(local_const)
        t_juncs = float(majiq_config.nrandom_junctions) / (50 * majiq_config.num_final_chunks)

        prb = min(1.0, float(t_juncs) / njuncs) * 100
        kk = np.random.choice(100, njuncs)
        indx = np.arange(njuncs)[kk <= prb]
        sample_junc = list(local_const)[0].get_coverage(ind_list)

        r_junctions = np.zeros(shape=(len(indx), sample_junc.shape[0], sample_junc.shape[1]))
        #r_junctions_gc = np.zeros(shape=(len(indx), list(local_const)[0].shape[0], list(local_const)[0].shape[1]))

        r_junctions_gc = []

        for jidx, jn in enumerate(np.array(list(local_const))[indx]):
            r_junctions[jidx, :, :] = jn.get_coverage(ind_list).toarray()
            if majiq_config.gcnorm:
                gc_array = jn.get_gc_content().toarray()
                gc_array = np.reshape(gc_array, (gc_array.shape[1],))
                r_junctions_gc.append([gc_vfunc[exp_idx](gc_array) for exp_idx in ind_list])
        if majiq_config.gcnorm:
            r_junctions = np.multiply(r_junctions, np.array(r_junctions_gc))

        qm = QueueMessage(QUEUE_MESSAGE_BUILD_CONST_JUNCTION, [r_junctions, name], chnk)
        out_queue.put(qm, block=True)

        for ss in dummy[name][0]:
            for st in dummy[name][1]:
                if ss.contained(st):
                    break
            else:
                qm = QueueMessage(QUEUE_MESSAGE_BUILD_LSV, [majiq_lsv.Queue_Lsv(ss, name), name], chnk)
                out_queue.put(qm, block=True)

        for st in dummy[name][1]:
            for ss in dummy[name][0]:
                if st.contained(ss):
                    break
            else:
                qm = QueueMessage(QUEUE_MESSAGE_BUILD_LSV, [majiq_lsv.Queue_Lsv(st, name), name], chnk)
                out_queue.put(qm, block=True)


