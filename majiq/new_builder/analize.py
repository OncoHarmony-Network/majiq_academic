import numpy as np

from majiq.grimoire.lsv import SSOURCE, STARGET, InvalidLSV
import majiq.src.config as majiq_config
import majiq.src.filter as majiq_filter
import majiq.grimoire.lsv as majiq_lsv

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


def lsv_detection(gn, only_real_data=False, out_queue=None, logging=None):

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
        r_junctions = np.array(list(local_const))[indx]

        for jix, jn in enumerate(r_junctions):
            out_queue.put([1, jn.get_coverage(ind_list), jn.get_gc_content(), name], block=True)

        for ss in dummy[name][0]:
            for st in dummy[name][1]:
                if ss.contained(st):
                    break
            else:
                out_queue.put([0, majiq_lsv.Queue_Lsv(ss, name), name], block=True)

        for st in dummy[name][1]:
            for ss in dummy[name][0]:
                if st.contained(ss):
                    break
            else:
                out_queue.put([0, majiq_lsv.Queue_Lsv(st, name), name], block=True)


