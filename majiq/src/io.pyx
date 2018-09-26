import datetime

import os
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector

from majiq.src.internals.grimoire cimport Gene, Exon, Junction
from majiq.src.internals import quant_lsv
from majiq.src.internals.psi cimport psi_distr_t, pair_int_t
from majiq.src.internals.qLSV cimport qLSV
from majiq.src.internals.io_utils cimport get_aggr_coverage

from majiq.src.gff import parse_gff3
from majiq.src.constants import *

from cython.parallel import prange

import pickle
import numpy as np
cimport numpy as np

# ctypedef vector[float] psi_distr_t ;

cdef list accepted_transcripts = ['mRNA', 'transcript', 'lnc_RNA', 'miRNA', 'ncRNA',
                                  'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'tRNA', 'pseudogenic_transcript']
cdef str transcript_id_keys = 'ID'
cdef list accepted_genes = ['gene', 'ncRNA_gene', 'pseudogene']
cdef list gene_name_keys = ['Name', 'gene_name']
cdef list gene_id_keys = ['ID', 'gene_id']


cdef int  read_gff(str filename, map[string, Gene*] all_genes, vector[string] gid_vec, object logging) except -1:
    """
    :param filename: GFF input filename
    :param list_of_genes: List of genes that will be updated with all the gene_id detected on the gff file
    :param logging: logger object
    :return: :raise RuntimeError:
    """

    cdef dict trcpt_id_dict = {}
    cdef dict exon_dict = {}

    cdef string chrom
    cdef char strand
    cdef int start, end
    cdef bint bb
    cdef list ind_list, tlist
    cdef string gene_id, key, gene_name, parent_tx_id
    # cdef map[string, Gene*] all_genes

    for record in parse_gff3(filename):
        chrom = record.seqid.encode('utf-8')

        strand = <char> record.strand.encode('UTF-8')[0]
        start = record.start
        end = record.end
        if record.type in accepted_genes:
            for gname_k in gene_name_keys:
                try:
                    gene_name = record.attributes[gname_k].encode('utf-8')
                    break
                except KeyError:
                    continue
            else:
                logging.info("Error, Gene doesn't contain one of the Name attribute  information values: "
                             "%s" % gene_name_keys)
            for gid_k in gene_id_keys:
                try:
                    gene_id = record.attributes[gid_k].encode('utf-8')
                    break
                except KeyError:
                    continue
            else:
                logging.info("Error, Gene doesn't contain one of the ID attribute information values: "
                             "%s" % gene_id_keys)
            if all_genes.count(gene_id)>0:
                raise RuntimeError('Two Different Genes with the same name %s' % gene_name)

            exon_dict[gene_id] = []
            all_genes[gene_id] = new Gene(gene_id, gene_name, chrom, strand, start, end)
            gid_vec.push_back(gene_id)
        elif record.type in accepted_transcripts:
            if transcript_id_keys not in record.attributes or 'Parent' not in record.attributes:
                logging.info("Error, Transcript doesn't contain one of the ID or parent attributes"
                             "information values: %s" % transcript_id_keys)
                continue
            transcript_name = record.attributes[transcript_id_keys]
            parent = record.attributes['Parent'].encode('utf-8')
            if all_genes.count(gene_id)==0:
                logging.error("Error, incorrect gff. mRNA %s doesn't have valid gene %s" % (transcript_name, parent))
                continue

            trcpt_id_dict[record.attributes['ID'].encode('utf-8')] = [parent, []]

        elif record.type == 'exon':
            parent_tx_id = record.attributes['Parent'].encode('utf-8')
            try:
                gene_id = trcpt_id_dict[parent_tx_id][0]
                exon_dict[gene_id].append((start, True))
                exon_dict[gene_id].append((end, False))
                trcpt_id_dict[parent_tx_id][1].append((start, end))

            except KeyError:
                logging.warning("Error, incorrect gff. exon "
                                "doesn't have valid mRNA %s" % parent_tx_id)

    for parent_tx_id, (gene_id, coord_list) in trcpt_id_dict.items():
        last_ss = FIRST_LAST_JUNC
        coord_list.sort(key=lambda x: (x[0], x[1]))
        # if gene_id == 'ENSMUSG00000006498': print (coord_list)
        for xx, yy in coord_list:
            key = ('%s-%s' % (last_ss, xx)).encode('utf-8')

            all_genes[gene_id].junc_map_[key] = new Junction(last_ss, xx, True)
            last_ss = yy

        key = ('%s-%s' % (last_ss, FIRST_LAST_JUNC)).encode('utf-8')
        all_genes[gene_id].junc_map_[key] = new Junction(last_ss, FIRST_LAST_JUNC, True)
    merge_exons(exon_dict, all_genes)
    return 0


cdef int merge_exons(dict exon_dict, map[string, Gene*]& all_genes) except -1:
    cdef list ex_list
    cdef string gne_id, key
    cdef tuple x
    cdef int ex_start, ex_end, nopen, coord
    cdef bint is_start


    for gne_id, ex_list in exon_dict.items():
        ex_list.sort(key=lambda x:(x[0], x[1]))
        ex_start = -1
        ex_end = -1
        nopen = 0
        for coord, is_start in ex_list:
            if is_start:
                if ex_end != -1:
                    key = ('%s-%s' % (ex_start, ex_end)).encode('utf-8')
                    all_genes[gne_id].exon_map_[key] = new Exon(ex_start, ex_end, True)
                    if nopen > 0 and (ex_end+4) < (coord-1):
                       pass
                        # all_genes[gne_id].create_annot_intron(ex_end+1, coord-1)
                        #tlist.append([ex_end+1, coord-1, 1, IR_TYPE])
                    ex_end = -1
                    nopen = 0
                    ex_start = coord

                ex_start = coord if ex_start == -1 or coord < ex_start else ex_start
                nopen += 1

            else:
                nopen -= 1
                ex_end = coord if coord > ex_end else ex_end

        if ex_end != -1:
            key = ('%s-%s' % (ex_start, ex_end)).encode('utf-8')
            all_genes[gne_id].exon_map_[key] = new Exon(ex_start, ex_end, True)


#######
# io API
#######
cdef int _load_db(str filename, object elem_dict, object genes_dict) except -1:
    cdef list names = ['id', 'name', 'chromosome', 'strand']

    with open(filename, 'rb') as fp:
        all_files = np.load(fp)
        genes_dict = {all_files['gene_info'][ii][0].decode('UTF-8'):{xx: all_files['gene_info'][ii][idx]
                                                                    for idx, xx in enumerate(names)}
                                                                    for ii in range(all_files['gene_info'].shape[0])}
        for xx in genes_dict.keys():
            elem_dict[xx] = all_files[xx]

cdef int dump_lsv_coverage(str filename, dict cov_dict, list type_list, list junc_info, str exp_name):
    dt=np.dtype('|S250, |S250')

    with open(filename, 'w+b') as ofp:
        cov_dict['lsv_types'] = np.array(type_list, dtype=dt)
        dt=np.dtype('|S250, u4, u4, f4, f4')
        cov_dict['junc_info'] = np.array(junc_info, dtype=dt)
        dt = np.dtype('|S250, |S25')
        cov_dict['meta'] = np.array([(exp_name, VERSION)], dtype=dt)
        np.savez(ofp, **cov_dict)


cdef int dump_lsv_coverage_mat(str filename, list cov_list, list type_list, list junc_info, str exp_name):
    dt=np.dtype('|S250, |S250')

    nlist = []
    xx = {}
    with open(filename, 'w+b') as ofp:

        xx['lsv_types'] = np.array(type_list, dtype=dt)
        dt=np.dtype('|S250, u4, u4, f4, f4')
        xx['junc_info'] = np.array(junc_info, dtype=dt)
        xx['coverage'] = np.array(cov_list, dtype=np.float32)
        dt = np.dtype('|S250, |S25')
        xx['meta'] = np.array([(exp_name, VERSION)], dtype=dt)
        np.savez(ofp, **xx)

cdef int _dump_elems_list(object elem_dict, object gene_info, str outDir) except -1:

    dt=np.dtype('|S250, |S250, |S32, S1, u4, u4')
    kk = [(xx['id'], xx['name'], xx['chromosome'], xx['strand'], xx['start'], xx['end']) for xx in gene_info.values()]
    elem_dict['gene_info'] = np.array(kk, dtype=dt)
    np.savez(get_build_temp_db_filename(outDir), **elem_dict)


cdef dict _get_extract_lsv_list(list list_of_lsv_id, list file_list):
    cdef dict result = {}
    cdef int n_exp = len(file_list)
    cdef str lsv_id, lsv_type, fname
    cdef int fidx
    cdef np.ndarray cov, lsv_cov

    for fidx, fname in enumerate(file_list):
        with open(fname, 'rb') as fp:
            data = np.load(fp)
            for lsv_id in list_of_lsv_id:
                try:
                    cov = data[lsv_id]
                except KeyError:
                    continue
                try:
                    # print(lsv_id, result[lsv_id].coverage[fidx].shape)
                    result[lsv_id].coverage[fidx] = cov
                except KeyError:
                    njunc = cov.shape[0]
                    msamples = cov.shape[1]
                    lsv_cov = np.zeros(shape=(n_exp, njunc, msamples),  dtype=float)
                    lsv_cov[fidx] = cov
                    result[lsv_id] = quant_lsv(lsv_id, lsv_cov)

    return result


cdef int dump_hettmp_file(str fname, np.ndarray[np.float32_t, ndim=2, mode="c"] osamps):
    with open(fname, 'w+b') as fp:
        np.save(fp, osamps)


cdef void get_coverage_mat_lsv(map[string, qLSV*]& result, list file_list, str weight_fname, int nthreads):
    cdef int n_exp = len(file_list)
    cdef str lid, lsv_type, fname
    cdef string lsv_id
    cdef int fidx, njunc, msamples, i
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] cov
    cdef dict weights
    cdef object data
    cdef int nlsv = result.size()
    cdef string empty_string = ''.encode('utf-8')
    cdef string prev_lsvid = empty_string
    cdef int prev_juncidx = 0

    # if weight_fname != "":
    #     weights = _load_weights(list_of_lsv_id, weight_fname)

    for fidx, fname in enumerate(file_list):
        print(fname)
        with open(fname, 'rb') as fp:
            fl = np.load(fp)
            data = fl['coverage']
            info = fl['junc_info']
            msamples = data.shape[1]
            idx = -1
            for row in info:
                idx += 1
                lsv_id = row[0]
                if result.count(lsv_id) > 0:
                    if prev_lsvid != lsv_id:
                        if prev_lsvid != empty_string:
                            result[prev_lsvid].add(<np.float32_t *> cov.data, msamples)
                        prev_lsvid = lsv_id
                        prev_juncidx = -1
                        njunc = result[lsv_id].get_num_ways()
                        cov = np.zeros(shape=(njunc, msamples), dtype=np.float32)
                    prev_juncidx += 1
                    cov[prev_juncidx] = data[idx]
            result[prev_lsvid].add(<np.float32_t *> cov.data, msamples)
    return



cdef void get_coverage_mat(map[string, vector[psi_distr_t]]& result, map[string, int] lsv_map, list file_list,
                            str weight_fname, int nthreads):
    cdef int n_exp = len(file_list)
    cdef str lid, lsv_type, fname
    cdef string lsv_id
    cdef int fidx, njunc, msamples, i
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] cov
    cdef dict weights
    cdef object data
    cdef int nlsv = lsv_map.size()
    cdef string empty_string = ''.encode('utf-8')
    cdef string prev_lsvid = empty_string
    cdef int prev_juncidx = 0

    # if weight_fname != "":
    #     weights = _load_weights(list_of_lsv_id, weight_fname)

    for fidx, fname in enumerate(file_list):
        print(fname)
        with open(fname, 'rb') as fp:
            fl = np.load(fp)
            data = fl['coverage']
            info = fl['junc_info']
            msamples = data.shape[1]
            idx = -1
            for row in info:
                idx += 1
                lsv_id = row[0]
                if lsv_map.count(lsv_id) > 0:
                    if prev_lsvid != lsv_id:
                        if prev_lsvid != empty_string:
                            # print (prev_lsvid)
                            get_aggr_coverage(result, prev_lsvid, <np.float32_t *> cov.data, njunc, msamples)

                        prev_lsvid = lsv_id
                        prev_juncidx = -1
                        njunc = lsv_map[lsv_id]
                        cov = np.zeros(shape=(njunc, msamples), dtype=np.float32)
                    prev_juncidx += 1
                    cov[prev_juncidx] = data[idx]
            get_aggr_coverage(result, prev_lsvid, <np.float32_t *> cov.data, njunc, msamples)
        print (fname, result.size())


cpdef void get_coverage_lsv(map[string, vector[psi_distr_t]]& result, vector[string] list_of_lsv_id, list file_list,
                            str weight_fname, int nthreads):
    # cdef map[string, vector[psi_distr_t]] result
    cdef int n_exp = len(file_list)
    cdef str lid, lsv_type, fname
    cdef string lsv_id
    cdef int fidx, njunc, msamples, i
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] cov
    cdef dict weights
    cdef object data
    cdef int nlsv = list_of_lsv_id.size()

    # if weight_fname != "":
    #     weights = _load_weights(list_of_lsv_id, weight_fname)

    for fidx, fname in enumerate(file_list):
        print(fname)
        with open(fname, 'rb') as fp:
            data = np.load(fp)
            for lsv_id in list_of_lsv_id:
                lid = lsv_id.decode('utf-8')
                try:
                    cov = data[lid]
                    njunc = cov.shape[0]
                    msamples = cov.shape[1]
                    with nogil:
                        get_aggr_coverage(result, lsv_id, <np.float32_t *> cov.data, njunc, msamples)

                except KeyError:
                        continue


cdef map[string, psi_distr_t] _load_weights(list lsv_list, str file_name):
    cdef dict out_dict = {}
    # cdef str file_name
    cdef str xx
    # file_name = get_weights_filename(outdir, name)
    with open(file_name, 'rb') as fp:
        all_wgts = np.load(fp)
        for xx in lsv_list:
            out_dict[xx] = all_wgts[xx]

    return out_dict

cdef list _extract_lsv_summary(list files, int minnonzero, int min_reads, dict types_dict, object junc_info,
                               list exp_name_list, dict epsi=None, int percent=-1, object logger=None):
    cdef dict lsv_types, lsv_list = {}
    cdef list lsv_id_list = []
    cdef int nfiles = len(files)
    cdef int fidx
    cdef str ff
    cdef dict lsv_junc_info = {}
    cdef np.ndarray mtrx, vals
    cdef np.ndarray jinfo

    if percent == -1:
        percent = nfiles / 2
        percent = percent + 1 if nfiles % 2 != 0 else percent
    percent = min(int(percent), nfiles)

    for fidx, ff in enumerate(files):
        if not os.path.isfile(ff):
            logger.error('File %s does not exist. Exiting execution' % ff)
            exit(-1)

        if logger:
            logger.info("Parsing file: %s" % ff)
        with open(ff, 'rb') as fp:
            all_files = np.load(fp)
            lsv_types = {yy[0]:[yy[1], 0] for yy in all_files['lsv_types']}
            jinfo = all_files['junc_info']
            xp = all_files['meta'][0]

            exp_name_list.append(xp[0])

            pre_lsv = jinfo[0][0]
            lsv_t = False
            epsi_t = []
            lsv_junc_info = {zz: [] for zz in lsv_types.keys()}

            for xx in jinfo:
                lsv_id = xx[0]
                lsv_junc_info[lsv_id].append([xx[1], xx[2]])
                lsv_types[lsv_id][1] += 1
                if xx[0] == pre_lsv:
                    lsv_t = lsv_t or (xx[3] >=min_reads and xx[4] >= minnonzero)
                    if epsi is not None:
                        epsi_t.append(xx[3])
                else:
                    pre_lsv = lsv_id
                    lsv_t = (xx[3] >=min_reads and xx[4] >= minnonzero)
                    try:
                        lsv_list[pre_lsv] += int(lsv_t)
                        if epsi is not None:
                            epsi[lsv_id] += np.array(epsi_t)
                    except KeyError:

                        lsv_list[pre_lsv] = int(lsv_t)
                        if epsi is not None:
                            epsi[lsv_id] = np.array(epsi_t)
                    epsi_t = []

        junc_info.update(lsv_junc_info)
        types_dict.update(lsv_types)

    if epsi is not None:
        for xx in epsi.keys():
            epsi[xx] = epsi[xx] / nfiles
            epsi[xx] = epsi[xx] / epsi[xx].sum()
            epsi[xx][np.isnan(epsi[xx])] = 1.0 / nfiles

    for xx, yy in lsv_list.items():
        if yy >= percent:
            lsv_id_list.append(xx)
        junc_info[xx] = np.array(junc_info[xx])

    return lsv_id_list


cdef list _extract_lsv_summary_old(list files, int minnonzero, int min_reads, object types_dict, object junc_info,
                               list exp_name_list, dict epsi=None, int percent=-1, object logger=None):
    cdef dict lsv_types, lsv_list = {}
    cdef list lsv_id_list = []
    cdef int nfiles = len(files)
    cdef int fidx
    cdef str ff
    cdef dict lsv_junc_info = {}
    cdef np.ndarray mtrx, vals
    cdef np.ndarray jinfo

    if percent == -1:
        percent = nfiles / 2
        percent = percent + 1 if nfiles % 2 != 0 else percent
    percent = min(int(percent), nfiles)

    for fidx, ff in enumerate(files):
        if not os.path.isfile(ff):
            logger.error('File %s does not exist. Exiting execution' % ff)
            exit(-1)

        if logger:
            logger.info("Parsing file: %s" % ff)
        with open(ff, 'rb') as fp:
            all_files = np.load(fp)
            lsv_types = {yy[0].decode('UTF-8'):yy[1].decode('UTF-8') for yy in all_files['lsv_types']}
            jinfo = all_files['junc_info']
            xp = all_files['meta'][0]

            exp_name_list.append(xp[0].decode('UTF-8'))

            pre_lsv = jinfo[0][0].decode('UTF-8')
            lsv_t = False
            epsi_t = []
            lsv_junc_info = {zz: [] for zz in lsv_types.keys()}

            for xx in jinfo:
                lsv_id = xx[0].decode('UTF-8')
                lsv_junc_info[lsv_id].append([xx[1], xx[2]])

                if xx[0] == pre_lsv:
                    lsv_t = lsv_t or (xx[3] >=min_reads and xx[4] >= minnonzero)
                    if epsi is not None:
                        epsi_t.append(xx[3])
                else:
                    pre_lsv = lsv_id
                    lsv_t = (xx[3] >=min_reads and xx[4] >= minnonzero)
                    try:
                        lsv_list[pre_lsv] += int(lsv_t)
                        if epsi is not None:
                            epsi[lsv_id] += np.array(epsi_t)
                    except KeyError:

                        lsv_list[pre_lsv] = int(lsv_t)
                        if epsi is not None:
                            epsi[lsv_id] = np.array(epsi_t)
                    epsi_t = []

        junc_info.update(lsv_junc_info)
        types_dict.update(lsv_types)

    if epsi is not None:
        for xx in epsi.keys():
            epsi[xx] = epsi[xx] / nfiles
            epsi[xx] = epsi[xx] / epsi[xx].sum()
            epsi[xx][np.isnan(epsi[xx])] = 1.0 / nfiles

    for xx, yy in lsv_list.items():
        if yy >= percent:
            lsv_id_list.append(xx)
        junc_info[xx] = np.array(junc_info[xx])

    return lsv_id_list




# cdef list _extract_lsv_summary_C(list files, int minnonzero, int min_reads,  map[string, string] types_dict,
#                                  object junc_info, list exp_name_list, dict epsi=None, int percent=-1, object logger=None):
#     cdef dict lsv_types, lsv_list = {}
#     cdef list lsv_id_list = []
#     cdef int nfiles = len(files)
#     cdef int fidx
#     cdef str ff
#     cdef dict lsv_junc_info = {}
#     cdef np.ndarray mtrx, vals
#     cdef np.ndarray jinfo
#
#     if percent == -1:
#         percent = nfiles / 2
#         percent = percent + 1 if nfiles % 2 != 0 else percent
#     percent = min(int(percent), nfiles)
#
#     for fidx, ff in enumerate(files):
#         if not os.path.isfile(ff):
#             logger.error('File %s does not exist. Exiting execution' % ff)
#             exit(-1)
#
#         if logger:
#             logger.info("Parsing file: %s" % ff)
#         with open(ff, 'rb') as fp:
#             all_files = np.load(fp)
#             lsv_types = {yy[0].decode('UTF-8'):yy[1].decode('UTF-8') for yy in all_files['lsv_types']}
#             jinfo = all_files['junc_info']
#             xp = all_files['meta'][0]
#
#             exp_name_list.append(xp[0].decode('UTF-8'))
#
#             pre_lsv = jinfo[0][0].decode('UTF-8')
#             lsv_t = False
#             epsi_t = []
#             lsv_junc_info = {zz: [] for zz in lsv_types.keys()}
#
#             for xx in jinfo:
#                 lsv_id = xx[0].decode('UTF-8')
#                 lsv_junc_info[lsv_id].append([xx[1], xx[2]])
#
#                 if xx[0] == pre_lsv:
#                     lsv_t = lsv_t or (xx[3] >=min_reads and xx[4] >= minnonzero)
#                     if epsi is not None:
#                         epsi_t.append(xx[3])
#                 else:
#                     pre_lsv = lsv_id
#                     lsv_t = (xx[3] >=min_reads and xx[4] >= minnonzero)
#                     try:
#                         lsv_list[pre_lsv] += int(lsv_t)
#                         if epsi is not None:
#                             epsi[lsv_id] += np.array(epsi_t)
#                     except KeyError:
#
#                         lsv_list[pre_lsv] = int(lsv_t)
#                         if epsi is not None:
#                             epsi[lsv_id] = np.array(epsi_t)
#                     epsi_t = []
#
#         junc_info.update(lsv_junc_info)
#         types_dict.update(lsv_types)
#
#     if epsi is not None:
#         for xx in epsi.keys():
#             epsi[xx] = epsi[xx] / nfiles
#             epsi[xx] = epsi[xx] / epsi[xx].sum()
#             epsi[xx][np.isnan(epsi[xx])] = 1.0 / nfiles
#
#     for xx, yy in lsv_list.items():
#         if yy >= percent:
#             lsv_id_list.append(xx)
#         junc_info[xx] = np.array(junc_info[xx])
#
#     return lsv_id_list

####
# API
##

cpdef tuple extract_lsv_summary(list files, int minnonzero, int min_reads, dict types_dict, dict junc_info,
                                dict epsi=None, int percent=-1, object logger=None):
    cdef list r
    cdef list exp_list = []
    r = _extract_lsv_summary(files, minnonzero, min_reads, types_dict, junc_info, exp_list, epsi, percent, logger)

    return r, exp_list


# cpdef  extract_lsv_summary_C(list files, int minnonzero, int min_reads, map[string, string] types_dict,
#                              map[string, np.ndarray] junc_info, dict epsi=None, int percent=-1, object logger=None):
#     cdef list r
#     cdef list exp_list = []
#     r = _extract_lsv_summary_C(files, minnonzero, min_reads, types_dict, junc_info, exp_list, epsi, percent, logger)
#
#     return r, exp_list


cpdef int add_elements_mtrx(dict new_elems, object shared_elem_dict):

    for gne, mtrx in new_elems.items():
        kk = shared_elem_dict[gne]
        shared_elem_dict[gne] = kk + new_elems[gne]



cpdef load_bin_file(filename, logger=None):
    if not os.path.exists(filename):
        if logger:
            logger.error('Path %s for loading does not exist' % filename)
        return

    fop = open(filename, 'rb')
    fast_pickler = pickle.Unpickler(fop)
    # fast_pickler.fast = 1
    data = fast_pickler.load()
    fop.close()
    return data


cpdef dump_bin_file(data, str filename):
    with open(filename, 'wb') as ofp:
        fast_pickler = pickle.Pickler(ofp, protocol=2)
        # fast_pickler.fast = 1
        fast_pickler.dump(data)

cpdef dict get_extract_lsv_list(list list_of_lsv_id, list file_list, bint aggr=True, str filename=None):

    if aggr:
        return  _get_extract_lsv_list(list_of_lsv_id, file_list)
    else:
        return _get_extract_lsv_list(list_of_lsv_id, file_list)


cpdef load_db(str filename, object elem_dict, object genes_dict):
    _load_db(filename, elem_dict, genes_dict)

cpdef dump_db(object genes_dict, object elem_dict, str outDir):
    _dump_elems_list(elem_dict, genes_dict, outDir)


cpdef store_weights(list lsv_list, np.ndarray wgts, str outdir, str name):
    file_name = get_weights_filename(outdir, name)
    dd = {xx:wgts[idx] for idx, xx in enumerate(lsv_list)}
    with open(file_name, 'w+b') as ofp:
        np.savez(ofp, **dd)

cpdef load_weights(list lsv_list, str outdir, str name):
    cdef dict out_dict = {}
    cdef str file_name
    cdef str xx

    file_name = get_weights_filename(outdir, name)
    with open(file_name, 'rb') as fp:
        all_wgts = np.load(fp)
        for xx in lsv_list:
            out_dict[xx] = all_wgts[xx]

    return out_dict
