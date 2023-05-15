import os
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector

from rna_majiq.src.internals.grimoire cimport Gene, Exon, Junction, coord_key_t
from rna_majiq.src.internals import quant_lsv
from rna_majiq.src.internals.qLSV cimport qLSV

from rna_majiq.src.gff import parse_gff3
import rna_majiq.src.constants as constants
from math import ceil
from cython.parallel import prange

import pickle
import numpy as np
cimport numpy as np



cdef list accepted_transcripts = ['mRNA', 'transcript', 'lnc_RNA', 'miRNA', 'ncRNA',
                                  'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'tRNA', 'pseudogenic_transcript',
                                  'C_gene_segment', 'D_gene_segment', 'J_gene_segment',
                                  'V_gene_segment', 'unconfirmed_transcript', 'three_prime_overlapping_ncrna']
cdef str transcript_id_keys = 'ID'
cdef list accepted_genes = ['gene', 'ncRNA_gene', 'pseudogene', 'ncRNA_gene', 'bidirectional_promoter_lncRNA']
cdef list gene_name_keys = ['Name', 'gene_name']
cdef list gene_id_keys = ['ID', 'gene_id']



cdef int  read_gff(str filename, map[string, Gene*] all_genes, vector[string] gid_vec, bint simpl, bint enable_anot_ir,
                   object logging, sqlite3 * db) except -1:
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
    cdef string gene_id, gene_name, parent_tx_id
    cdef coord_key_t key
    # cdef map[string, Gene*] all_genes

    for record in parse_gff3(filename):

        if record.strand is None or record.seqid is None:
            continue

        chrom = record.seqid.encode('utf-8')

        strand = <char> record.strand.encode('UTF-8')[0]
        start = record.start
        end = record.end
        if record.type in accepted_genes:
            for gid_k in gene_id_keys:
                try:
                    gene_id = record.attributes[gid_k].encode('utf-8')
                    break
                except KeyError:
                    continue
            else:
                logging.warning(
                    "Error, gene (record.type=%s) doesn't have an attribute"
                    " recognized for inferring gene_id (%s)."
                    " (record.attributes=%s)",
                    record.type,
                    gene_id_keys,
                    record.attributes,
                )
                continue  # we cannot process this gene without the gene_id
            if all_genes.count(gene_id)>0:
                raise RuntimeError(
                    "Two gene GFF3 records with the same gene_id %s", gene_id
                )

            for gname_k in gene_name_keys:
                try:
                    gene_name = record.attributes[gname_k].encode('utf-8')
                    break
                except KeyError:
                    continue
            else:
                gene_name = gene_id  # gene_id is fallback for no gene_name

            exon_dict[gene_id] = []
            all_genes[gene_id] = new Gene(gene_id, gene_name, chrom, strand, start, end)
            gid_vec.push_back(gene_id)
        elif record.type in accepted_transcripts:
            if transcript_id_keys not in record.attributes or 'Parent' not in record.attributes:
                logging.info(
                    "Error, transcript (record.type=%s) doesn't have an"
                    " attribute recognized for inferring transcript_id (%s)"
                    " or transcript parent (Parent). (record.attributes=%s)",
                    record.type,
                    transcript_id_keys,
                    record.attributes,
                )
                continue
            transcript_name = record.attributes[transcript_id_keys]
            parent = record.attributes['Parent'].encode('utf-8')
            if all_genes.count(gene_id)==0:
                logging.error(
                    "Error, incorrect gff. transcript %s doesn't have valid gene %s",
                    transcript_name,
                    parent,
                )
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
                logging.warning(
                    "Error, incorrect gff. exon doesn't have valid transcript %s",
                    parent_tx_id,
                )
    # end loop over records in GFF3 file

    for parent_tx_id, (gene_id, coord_list) in trcpt_id_dict.items():

        last_ss = constants.FIRST_LAST_JUNC
        coord_list.sort(key=lambda x: (x[0], x[1]))
        if len(coord_list) == 0 : continue
        for xx, yy in coord_list:

            sg_transcript_exon(db, gene_id, parent_tx_id, xx, yy)

            key = coord_key_t(last_ss, xx)

            if all_genes[gene_id].junc_map_.count(key) == 0:
                all_genes[gene_id].junc_map_[key] = new Junction(last_ss, xx, True, simpl)
            last_ss = yy

        key = coord_key_t(last_ss, constants.FIRST_LAST_JUNC)
        if all_genes[gene_id].junc_map_.count(key) == 0:
            all_genes[gene_id].junc_map_[key] = new Junction(last_ss, constants.FIRST_LAST_JUNC, True, simpl)
    merge_exons(exon_dict, all_genes, simpl, enable_anot_ir)
    return 0


cdef int merge_exons(dict exon_dict, map[string, Gene*]& all_genes, bint simpl, bint enable_anot_ir) except -1:
    cdef list ex_list
    cdef string gne_id
    cdef coord_key_t key
    cdef tuple x
    cdef int ex_start, ex_end, nopen, coord
    cdef bint is_start


    for gne_id, ex_list in exon_dict.items():
        ex_list.sort(key=lambda x:(x[0], not x[1]))
        ex_start = constants.EMPTY_COORD
        ex_end = constants.EMPTY_COORD
        nopen = 0

        for coord, is_start in ex_list:
            if is_start:
                if ex_end != constants.EMPTY_COORD:
                    start1 = ex_end -10  if ex_start == constants.EMPTY_COORD else ex_start
                    end1 = ex_start +10  if ex_end == constants.EMPTY_COORD else ex_end
                    key = coord_key_t(start1, end1)
                    all_genes[gne_id].exon_map_[key] = new Exon(ex_start, ex_end, True)

                    if nopen > 0 and ex_end < coord:
                        # create annotated intron (function adjusts coordinates using exon coordinates)
                        all_genes[gne_id].create_annot_intron(ex_end, coord, simpl, enable_anot_ir)
                    ex_end = constants.EMPTY_COORD
                    ex_start = coord

                ex_start = coord if ex_start == constants.EMPTY_COORD or coord < ex_start else ex_start
                nopen += 1

            else:
                nopen -= 1
                ex_end = coord if coord > ex_end else ex_end

        if ex_end != constants.EMPTY_COORD:
            key = coord_key_t(ex_start, ex_end)
            all_genes[gne_id].exon_map_[key] = new Exon(ex_start, ex_end, True)


#######
# io API
#######

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
        xx['meta'] = np.array([(exp_name, constants.VERSION)], dtype=dt)
        np.savez(ofp, **xx)


cdef int dump_hettmp_file(str fname, np.ndarray[np.float32_t, ndim=2, mode="c"] osamps):
    with open(fname, 'w+b') as fp:
        np.save(fp, osamps)


cdef void get_coverage_mat_lsv(map[string, qLSV*]& result, list file_list, int nthreads,
                               bint fltr, int minreads, int minnonzero):
    """ Add coverage from file_list majiq files to named LSVs in result

    Parameters
    ----------
    result: map[string, qLSV*]
        Map from lsv_id string to LSV object pointer with coverage information
        and enabled information, updated by reference
    file_list: List[str]
        List of file paths with MAJIQ coverage information
    nthreads: int
        Number of threads, which is currently unused
    fltr: bool
        If True, use all LSVs regardless of coverage. If False (used by HET),
        only add coverage for an experiment if passes minreads/minnonzero
        filters.
    minreads, minnonzero: int
        per-experiment filters to be used on reads/positions

    Notes
    -----
    Each LSV has a global flag indicating whether it is enabled, which is
    modified by qLSV::set_bool(). This is used by HET on a per-experiment
    basis. If fltr is True, any LSV found in any of the MAJIQ files that also
    passed group filters (and is thus in result) will have this set to enabled.
    Otherwise (fltr is False), if the LSV was previously enabled, it will stay
    enabled, but otherwise, it will only be set to enabled if it passes
    per-experiment filters in at least one of the input MAJIQ files. This is
    used in HET to disable LSVs on per-experiment basis by disabling all LSVs
    and only re-enabling the LSVs found above per-experiment filters in each
    file individually each time (file_list of length 1 each time after
    resetting all LSVs to disabled)
    """
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
    cdef bint bflt = not fltr

    for fidx, fname in enumerate(file_list):
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
                            # update add coverage/set status for previous LSV
                            result[prev_lsvid].add(<np.float32_t *> cov.data, msamples)
                            result[prev_lsvid].set_bool(bflt or result[prev_lsvid].is_enabled())
                        # if fltr, LSV will be used no matter what
                        bflt = fltr or (row[3] >=minreads and row[4] >= minnonzero)
                        prev_lsvid = lsv_id
                        prev_juncidx = -1


                        njunc = result[lsv_id].get_num_ways()
                        cov = np.zeros(shape=(njunc, msamples), dtype=np.float32)
                    else:
                        bflt = bflt or (row[3] >=minreads and row[4] >= minnonzero)
                    prev_juncidx += 1
                    cov[prev_juncidx] = data[idx]
            # make sure to add coverage/set status for final LSV
            result[prev_lsvid].add(<np.float32_t *> cov.data, msamples)
            result[prev_lsvid].set_bool(bflt or result[prev_lsvid].is_enabled())

    return


cdef list _extract_lsv_summary(list files, int minnonzero, int min_reads, dict types_dict, object junc_info,
                               list exp_name_list, dict o_epsi=None, dict prior_conf=None, float nexp=-1,
                               object logger=None):
    cdef dict lsv_types, lsv_list = {}, lsv_list_prior = {}
    cdef list lsv_id_list = []
    cdef int nfiles = len(files)
    cdef int fidx
    cdef str ff
    cdef dict lsv_junc_info = {}
    cdef np.ndarray mtrx, vals
    cdef np.ndarray jinfo
    cdef dict epsi = {}
    cdef int percent

    percent = min(nfiles, int(ceil(nfiles * nexp if nexp < 1 else nexp)))

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
            lsv_t_prior = False
            epsi_t = []
            lsv_junc_info = {zz: [] for zz in lsv_types.keys()}

            for xx in jinfo:
                lsv_id = xx[0]

                lsv_junc_info[lsv_id].append([xx[1], xx[2]])
                lsv_types[lsv_id][1] += 1

                if lsv_id == pre_lsv:
                    lsv_t = lsv_t or (xx[3] >=min_reads and xx[4] >= minnonzero)
                    if prior_conf is not None:
                        lsv_t_prior = lsv_t_prior or (xx[3] >=prior_conf['mreads'] and xx[4] >= prior_conf['mpos'])
                    epsi_t.append(xx[3])
                else:
                    try:
                        lsv_list[pre_lsv] += int(lsv_t)
                        if o_epsi is not None:
                            lsv_list_prior[pre_lsv] += int(lsv_t_prior)
                            epsi[pre_lsv] += np.array(epsi_t)

                    except KeyError:
                        lsv_list[pre_lsv] = int(lsv_t)
                        if o_epsi is not None:
                            lsv_list_prior[pre_lsv] = int(lsv_t_prior)
                            epsi[pre_lsv] = np.array(epsi_t)

                    epsi_t = [xx[3]]
                    pre_lsv = lsv_id
                    lsv_t = (xx[3] >=min_reads and xx[4] >= minnonzero)
                    if prior_conf is not None:
                        lsv_t_prior = lsv_t or (xx[3] >=prior_conf['mreads'] and xx[4] >= prior_conf['mpos'])
            try:
                lsv_list[pre_lsv] += int(lsv_t)
                if o_epsi is not None:
                    lsv_list_prior[pre_lsv] += int(lsv_t_prior)
                    epsi[pre_lsv] += np.array(epsi_t)
            except KeyError:
                lsv_list[pre_lsv] = int(lsv_t)
                if o_epsi is not None:
                    lsv_list_prior[pre_lsv] = int(lsv_t_prior)
                    epsi[pre_lsv] = np.array(epsi_t)

        junc_info.update(lsv_junc_info)
        types_dict.update(lsv_types)

    if o_epsi is not None:
        for xx in epsi.keys():
            if lsv_list_prior[xx] >= percent:
                o_epsi[xx] = epsi[xx] / nfiles
                o_epsi[xx] = o_epsi[xx] / o_epsi[xx].sum()
                o_epsi[xx][np.isnan(o_epsi[xx])] = 1.0 / nfiles

    for xx, yy in lsv_list.items():
        if yy >= percent:
            lsv_id_list.append(xx)
        junc_info[xx] = np.array(junc_info[xx])

    return lsv_id_list


####
# API
##

cpdef tuple extract_lsv_summary(list files, int minnonzero, int min_reads, dict types_dict, dict junc_info,
                                dict epsi=None, dict prior_conf=None, float percent=-1, object logger=None):
    cdef list r
    cdef list exp_list = []
    r = _extract_lsv_summary(files, minnonzero, min_reads, types_dict, junc_info, exp_list, epsi, prior_conf,
                             percent, logger)

    return r, exp_list

