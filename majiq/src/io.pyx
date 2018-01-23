import datetime

# import h5py
import os
from majiq.grimoire.junction cimport Junction
from majiq.grimoire.junction import Junction
from majiq.grimoire.exon cimport Exon, Intron
from majiq.grimoire.exon import Exon
from majiq.grimoire.lsv import quant_lsv

from majiq.src.gff import parse_gff3
from majiq.src.constants import *

import pickle
import numpy as np
cimport numpy as np


cdef list accepted_transcripts = ['mRNA', 'transcript', 'lnc_RNA', 'miRNA', 'ncRNA',
                                  'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'tRNA', 'pseudogenic_transcript']
cdef str transcript_id_keys = 'ID'
cdef list accepted_genes = ['gene', 'ncRNA_gene', 'pseudogene']
cdef list gene_name_keys = ['Name', 'gene_name']
cdef list gene_id_keys = ['ID', 'gene_id']


cdef int _read_gff(str filename, str outDir, object elem_dict,  object all_genes, object logging=None) except -1:
    """
    :param filename: GFF input filename
    :param list_of_genes: List of genes that will be updated with all the gene_id detected on the gff file
    :param logging: logger object
    :return: :raise RuntimeError:
    """
    # cdef list all_genes = []
    cdef dict gene_id_dict = {}
    cdef dict trcpt_id_dict = {}
    cdef dict exon_dict = {}
    # cdef dict elem_dict = {}

    cdef str chrom, name, strand
    cdef int start, end
    cdef bint bb
    cdef list ind_list

    for record in parse_gff3(filename):
        chrom = record.seqid
        strand = record.strand
        start = record.start
        end = record.end

        if record.type in accepted_genes:
            for gname_k in gene_name_keys:
                try:
                    gene_name = record.attributes[gname_k]
                    break
                except KeyError:
                    continue
            else:
                logging.info("Error, Gene doesn't contain one of the Name attribute  information values: "
                             "%s" % gene_name_keys)

            for gid_k in gene_id_keys:
                try:
                    gene_id = record.attributes[gid_k]
                    break
                except KeyError:
                    continue
            else:
                logging.info("Error, Gene doesn't contain one of the ID attribute information values: "
                             "%s" % gene_id_keys)

            exon_dict[gene_id] = []
            #all_genes.append((gene_id, gene_name, chrom, strand, start, end))
            all_genes[gene_id] = {'id':gene_id, 'name':gene_name, 'chromosome': chrom,
                                  'strand': strand, 'start': start, 'end':end, 'nreads': 0}

            #all_genes[chrom].append((gene_id, gene_name, chrom, strand, start, end))
            elem_dict[gene_id] = []

            if gene_id in gene_id_dict:
                raise RuntimeError('Two Different Genes with the same name %s' % gene_name)
            gene_id_dict[gene_id] = gene_id

        elif record.type in accepted_transcripts:
            if transcript_id_keys not in record.attributes or 'Parent' not in record.attributes:
                logging.info("Error, Transcript doesn't contain one of the ID or parent attributes"
                             "information values: %s" % transcript_id_keys)
                continue
            transcript_name = record.attributes[transcript_id_keys]
            parent = record.attributes['Parent']
            try:
                gn_id = gene_id_dict[parent]
                trcpt_id_dict[record.attributes['ID']] = [gn_id, []]

            except KeyError:
                logging.error("Error, incorrect gff. mRNA %s doesn't have valid gene %s" % (transcript_name, parent))
                #raise

        elif record.type == 'exon':
            parent_tx_id = record.attributes['Parent']
            try:
                gn_id = trcpt_id_dict[parent_tx_id][0]
                exon_dict[gn_id].append((start, True))
                exon_dict[gn_id].append((end, False))
                trcpt_id_dict[parent_tx_id][1].append((start, end))

            except KeyError:
                logging.warning("Error, incorrect gff. exon "
                                "doesn't have valid mRNA %s" % parent_tx_id)


    for parent_tx_id, (gn_id, coord_list) in trcpt_id_dict.items():
        last_ss = FIRST_LAST_JUNC
        coord_list.sort(key=lambda x: (x[0], x[1]))
        for xx, yy in coord_list:
            tlist = elem_dict[gn_id]
            #elem_dict[gn_id].append([last_ss, xx, 1, J_TYPE])
            tlist.append([last_ss, xx, 1, J_TYPE])
            elem_dict[gn_id] = tlist
            last_ss = yy

        elem_dict[gn_id].append([last_ss, FIRST_LAST_JUNC, 1, J_TYPE])
    merge_exons(exon_dict, elem_dict)


cdef int merge_exons(dict exon_dict, object elem_dict) except -1:
    cdef list ex_list
    cdef str gne_id
    cdef tuple x
    cdef int ex_start, ex_end, nopen, coord
    cdef bint is_start


    for gne_id, ex_list in exon_dict.items():
        ex_list.sort(key=lambda x:(x[0],x[1]))
        ex_start = -1
        ex_end = -1
        nopen = 0
        for coord, is_start in ex_list:
            if is_start:
                if ex_end != -1:
                    elem_dict[gne_id].append([ex_start, ex_end, 1, EX_TYPE])
                    if nopen > 0 and (ex_end+4) < (coord-1):
                        elem_dict[gne_id].append([ex_end+1, coord-1, 1, IR_TYPE])
                    ex_end = -1
                    nopen = 0
                    ex_start = coord

                ex_start = coord if ex_start == -1 or coord < ex_start else ex_start
                nopen += 1

            else:
                nopen -= 1
                ex_end = coord if coord > ex_end else ex_end

        if ex_end != -1:
            elem_dict[gne_id].append([ex_start, ex_end, 1, EX_TYPE])

        # _dump_elems_list(db_f, gne_id, np.array(elem_mtrx))
        #elem_dict[gne_id] = np.array(elem_dict[gne_id])


#######
# HDF5 API
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

cdef int _dump_lsv_coverage(str filename, dict cov_dict, list type_list):
    dt=np.dtype('|S250, |S250')
    with open(filename, 'w+b') as ofp:
        cov_dict['lsv_types'] = np.array(type_list, dtype=dt)
        np.savez(ofp, **cov_dict)

cdef int _dump_elems_list(object elem_dict, object gene_info, str outDir) except -1:

    dt=np.dtype('|S250, |S250, |S32, S1, u4, u4')
    kk = [(xx['id'], xx['name'], xx['chromosome'], xx['strand'], xx['start'], xx['end']) for xx in gene_info.values()]
    elem_dict['gene_info'] = np.array(kk, dtype=dt)
    np.savez(get_build_temp_db_filename(outDir), **elem_dict)


cdef int _read_junction(list row, str gne_id, dict jjs, list exs, list irs, int default_index) except -1:
    jjs[(row[0], row[1])] = Junction(row[0], row[1], gne_id, default_index, annot=bool(row[2]))


cdef int _read_exon(list row, str gne_id, dict jjs, list exs, list irs, int default_index) except -1:
    exs.append(Exon(row[0], row[1], annot=bool(row[2])))


cdef int _read_ir(list row, str gne_id, dict jjs, list exs, list irs, int default_index) except -1:
    irs.append(Intron(row[0], row[1], annot=bool(row[2]), db_idx=-1))


cdef int _pass_ir(list row, str gne_id, dict jjs, list exs, list irs, int default_index) except -1:
    pass


cdef _get_extract_lsv_list(list list_of_lsv_id, list file_list):
    cdef list result = []
    cdef int n_exp = len(file_list)
    cdef str lsv_id, lsv_type, fname
    cdef int fidx
    cdef np.ndarray cov, lsv_cov

    for lsv_id in list_of_lsv_id:
        lsv_type = None
        lsv_cov = None
        for fidx, fname in enumerate(file_list):
            with open(fname, 'rb') as fp:
                try:
                    cov = np.load(fp)[lsv_id][:,:-2]
                except KeyError:
                    continue
            njunc = cov.shape[0]
            msamples = cov.shape[1]

            if lsv_cov is None:
                lsv_cov = np.zeros(shape=(n_exp, njunc, msamples),  dtype=float)

            lsv_cov[fidx] = cov
        qq = quant_lsv(lsv_id, lsv_cov)
        result.append(qq)
    return result


####
# API
##

# def read_meta_info(list_of_files):
#     meta = {'experiments': []}
#     m_samples = None
#     for fl in list_of_files:
#         with h5py.File(fl, 'r') as fp :
#             if m_samples is not None:
#                 assert m_samples == fp.attrs['m_samples'], "uneven number of bootstrap samples"
#             else:
#                 m_samples = fp.attrs['m_samples']
#             meta['experiments'].append(fp.attrs['sample_id'])
#             try:
#                 if meta['genome'] != fp.attrs['genome']:
#                     raise RuntimeError('Combining experiments from different genome assemblies. Exiting')
#             except KeyError:
#                 meta['genome'] = fp.attrs['genome']
#                 continue
#     meta['m_samples'] = m_samples
#     return meta

cpdef extract_lsv_summary(list files, int minnonzero, int min_reads, object types_dict, dict epsi=None,
                          int percent=-1, object logger=None):
    cdef dict lsv_list = {}
    cdef dict lsv_types = {}
    cdef int nfiles = len(files)
    cdef int fidx
    cdef str ff, xx
    cdef np.ndarray mtrx, vals

    if percent == -1:
        percent = nfiles / 2
        percent = percent + 1 if nfiles % 2 != 0 else percent
    percent = min(int(percent), nfiles)

    for fidx, ff in enumerate(files):
        if logger:
            logger.info("Parsing file: %s" % ff)

        with open(ff, 'rb') as fp:
            all_files = np.load(fp)
            # print ([xx for xx in all_files['lsv_types']])
            lsv_types = {all_files['lsv_types'][ii][0].decode('UTF-8'):all_files['lsv_types'][ii][1].decode('UTF-8')
                         for ii in range(all_files['lsv_types'].shape[0])}
            for xx in lsv_types:
                lsv_data = all_files[xx][:,-2:]

                try:
                    lsv_list[xx] += int(np.any(np.logical_and(lsv_data[:, 0] >= minnonzero, lsv_data[:,1] >= min_reads)))
                    if epsi is not None:
                        epsi[xx] += lsv_data[:, 0]

                except KeyError:
                    lsv_list[xx] = int(np.any(np.logical_and(lsv_data[:, 0] >= minnonzero, lsv_data[:,1] >= min_reads)))
                    if epsi is not None:
                        epsi[xx] = lsv_data[:, 0]

        types_dict.update(lsv_types)

    if epsi is not None:
        for xx in epsi.keys():
            epsi[xx] = epsi[xx] / nfiles
            epsi[xx] = epsi[xx] / epsi[xx].sum()
            epsi[xx][np.isnan(epsi[xx])] = 1.0 / nfiles

    lsv_id_list = [xx for xx, yy in lsv_list.items() if yy >= percent]

    return lsv_id_list


cpdef int dump_lsv_coverage(out_f, cov_list, attrs_list):
    _dump_lsv_coverage(out_f, cov_list, attrs_list)

cpdef int parse_annot(str filename, str out_dir,  object elem_dict,  object all_genes, object logging=None):

    try:
        _read_gff(filename=filename, outDir=out_dir, elem_dict=elem_dict, all_genes=all_genes, logging=logging)
    except Exception as e:
        print e
        raise


cpdef int from_matrix_to_objects( str gne_id, object elem_dicts, dict dict_junctions,
                                 list list_exons, list list_introns=None, int default_index=-1):
    cdef dict func_list
    cdef list elem

    if list_introns is not None:
        func_list = {EX_TYPE: _read_exon, IR_TYPE: _read_ir, J_TYPE: _read_junction}
    else:
        func_list = {EX_TYPE: _read_exon, IR_TYPE: _pass_ir, J_TYPE: _read_junction}

    for elem in elem_dicts:
        func_list[elem[3]](elem, gne_id, dict_junctions, list_exons,
                             list_introns, default_index)


cpdef int add_elements_mtrx(dict new_elems, object shared_elem_dict):

    for gne, mtrx in new_elems.items():
        kk = shared_elem_dict[gne]
        shared_elem_dict[gne] = kk + new_elems[gne]


cpdef dump_elements(object genes_dict, object elem_dict, str outDir):
    _dump_elems_list(elem_dict, genes_dict, outDir)


def load_bin_file(filename, logger=None):
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


def dump_bin_file(data, str filename):
    with open(filename, 'wb') as ofp:
        fast_pickler = pickle.Pickler(ofp, protocol=2)
        # fast_pickler.fast = 1
        fast_pickler.dump(data)

def get_extract_lsv_list(list list_of_lsv_id, list file_list):
    return _get_extract_lsv_list(list_of_lsv_id, file_list)


def load_db(str filename, object elem_dict, object genes_dict):
    _load_db(filename, elem_dict, genes_dict)
