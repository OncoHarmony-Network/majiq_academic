"""
Functions to filter junction pairs by number of positions covered or number of reads
"""
from majiq.src.constants import *


def reliable_in_data(junc, exp_idx, minnonzero=2, min_reads=3):
    min_read_x_exp = min_reads
    min_npos_x_exp = minnonzero
    in_data_filter = False
    if not junc.all_data:
        numpos = np.count_nonzero(junc.coverage[exp_idx])
        totreads = junc.get_read_num(exp_idx)
    else:
        kk = junc.get_coverage()
        totreads = kk[0][exp_idx]
        numpos = kk[1][exp_idx]
        # numpos = np.count_nonzero(cov)
        # totreads = cov.sum()
    if totreads >= min_read_x_exp and numpos >= min_npos_x_exp:
        in_data_filter = True
    return in_data_filter


def filter_message(when, value, logger, junc):
    message = "%s (Filter=%s). %s" % (when, value, len(junc[0]))
    if logger:
        if type(logger) == bool:
            print(message)
        else:
            logger.debug(message)


def merge_files_hdf5(lsv_dict, lsv_summarized, minnonzero, min_reads, percent=-1, logger=None):

    logger.debug("Quantifible filter...")
    nfiles = lsv_summarized.shape[0]
    if percent == -1:
        percent = nfiles / 2
        percent = percent + 1 if nfiles % 2 != 0 else percent

    percent = min(percent, nfiles)
    vals = ((lsv_summarized[:, :, 0] >= minnonzero) * (lsv_summarized[:, :, 1] >= min_reads)).sum(axis=0) >= percent
    list_of_lsvs = [kk for kk, vv in lsv_dict.items() if vals[vv].sum() >= 1]
    return list_of_lsvs


def quantifiable_in_group_to_hdf5(hdf5_p, list_of_experiments, minnonzero, min_reads, effective_readlen,
                                  filter_vals=None, logger=None):
    logger.debug("Quantifible filter...")
    nexp = len(list_of_experiments)

    filtr = nexp / 2
    if nexp % 2 != 0:
        filtr += 1

    filt_exp = {}
    for idx, exp in enumerate(list_of_experiments):
        temp = lsv_quantifiable(exp, minnonzero, min_reads, filter_vals, logger)
        for ldx, lsv in enumerate(temp[1]):
            if not lsv[1] in filt_exp:
                filt_exp[lsv[0]] = 0
            filt_exp[lsv[0]] += 1

    tlb = {}
    info_tlb = {}
    info_tlb2 = {}
    for idx, exp in enumerate(list_of_experiments):
        for idx_lsv, lsv in enumerate(exp[1]):
            if not lsv[0] in tlb:
                tlb[lsv[0]] = [None] * nexp
                nways = len(exp[0][idx_lsv])
                info_tlb[lsv[0]] = lsv
                info_tlb2[lsv[0]] = np.zeros(shape=exp[0][idx_lsv].shape)
            tlb[lsv[0]][idx] = idx_lsv

    # filtered = []
    # filtered_info = []
    info = tlb.keys()
    old_shape = len(info)
    lsv_idx = 0

    hdf5_p.create_dataset(JUNCTIONS_DATASET_NAME, (len(info)*2, nexp, effective_readlen),
                          maxshape=(None, nexp, effective_readlen))

    list_lsv = []
    for ii in info:
        if ii not in filt_exp or filt_exp[ii] < filtr:
            continue

        type_lsv = info_tlb[ii][1]
        list_lsv.append(ii)
        lsv_idx += 1
        h_lsv = hdf5_p.create_group(ii)
        h_lsv.attrs['id'] = ii
        h_lsv.attrs['type'] = type_lsv
        info_tlb[ii][2].copy(info_tlb[ii][2], h_lsv, name='visuals')
        all_vals = []
        for idx, exp in enumerate(list_of_experiments):

            if tlb[ii][idx] is not None:
                val = exp[0][tlb[ii][idx]]
            else:
                val = info_tlb2[ii]
            all_vals.append(val)
        njunc = len(val)
        if lsv_idx + njunc > old_shape:
            shp = hdf5_p[JUNCTIONS_DATASET_NAME].shape
            shp_new = shp[0] + len(info)
            hdf5_p[JUNCTIONS_DATASET_NAME].resize((shp_new, shp[1]))

        hdf5_p[JUNCTIONS_DATASET_NAME][lsv_idx:lsv_idx+njunc] = np.reshape(np.array(all_vals), (njunc, nexp, effective_readlen))
        h_lsv.attrs['coverage'] = hdf5_p[JUNCTIONS_DATASET_NAME].regionref[lsv_idx:lsv_idx + njunc]

    return list_lsv


def quantifiable_in_group(list_of_experiments, minnonzero, min_reads, filter_vals=None, logger=None):
    nexp = len(list_of_experiments)

    filtr = nexp / 2
    if nexp % 2 != 0:
        filtr += 1

    filt_exp = {}
    for idx, exp in enumerate(list_of_experiments):
        temp = lsv_quantifiable(exp, minnonzero, min_reads, filter_vals, logger)
        for ldx, lsv in enumerate(temp[1]):
            if not lsv[1] in filt_exp:
                filt_exp[lsv[1]] = 0
            filt_exp[lsv[1]] += 1

    tlb = {}
    info_tlb = {}
    info_tlb2 = {}
    for idx, exp in enumerate(list_of_experiments):
        for idx_lsv, lsv in enumerate(exp[1]):
            if not lsv[1] in tlb:
                tlb[lsv[1]] = [None] * nexp
                nways = len(exp[0][idx_lsv])
                info_tlb[lsv[1]] = lsv
                info_tlb2[lsv[1]] = np.zeros(shape=exp[0][idx_lsv].shape)
            tlb[lsv[1]][idx] = idx_lsv

    filtered = []
    filtered_info = []
    info = tlb.keys()
    for ii in info:
        if not ii in filt_exp:
            continue
        pres = filt_exp[ii]
        if pres < filtr:
            continue
        lsv = []
        idlsv = info_tlb[ii]

        for idx, exp in enumerate(list_of_experiments):

            local_indx = tlb[ii][idx]
            if not local_indx is None:
                val = exp[0][local_indx]
            else:
                val = info_tlb2[ii]

            lsv.append(val)
        filtered.append(lsv)
        filtered_info.append(idlsv)

    return filtered, filtered_info


def lsv_quantifiable(list_lsv_tuple, minnonzero, min_reads, filter_vals=None, logger=False, fon=[True, True],
                     const=False):
    filter_message("Before quantifiable_filter", minnonzero, logger, np.array(list_lsv_tuple))
    filtered = []
    filtered_info = []
    k = 2
    if filter_vals is None:
        for lsvdx, lsv in enumerate(list_lsv_tuple[0]):
            for idx in range(lsv.shape[0]):
                if ((not fon[1] or np.count_nonzero(lsv[idx]) >= minnonzero) and
                        (not fon[0] or lsv[idx].sum() >= min_reads)):
                    filtered.append(list_lsv_tuple[0][lsvdx])
                    if not const:
                        filtered_info.append(list_lsv_tuple[1][lsvdx])
                    break
    else:
        for lsvdx, lsv in enumerate(filter_vals):
            if lsv:
                filtered.append(list_lsv_tuple[0][lsvdx])
                filtered_info.append(list_lsv_tuple[1][lsvdx])

    # filter_message("After quantifiable_filter", minnonzero, logger, array(filtered))
    return filtered, filtered_info


def lsv_intersection(lsv_list1, lsv_list2, bycol=False):
    lsv_match = [[], []]
    match_info = []

    ids1 = set([xx[1] for xx in lsv_list1[1]])
    ids2 = set([xx[1] for xx in lsv_list2[1]])
    matched_names = ids1.intersection(ids2)

    for ii in matched_names:
        typ_check = ''
        for idx, nm in enumerate(lsv_list1[1]):
            if nm[1] == ii:
                typ_check = nm[2]
                lsv_match[0].append(lsv_list1[0][idx])
                match_info.append(nm)
                break
        for idx, nm in enumerate(lsv_list2[1]):

            if nm[1] == ii:
                if typ_check != nm[2]:
                    continue
                lsv_match[1].append(lsv_list2[0][idx])
                break

    return lsv_match, match_info

