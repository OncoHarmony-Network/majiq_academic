import datetime
cimport numpy as np
import numpy as np
from numpy cimport ndarray

import h5py
import os
from majiq.grimoire.junction cimport Junction
from majiq.grimoire.junction import Junction
from majiq.grimoire.exon cimport Exon, Intron
from majiq.grimoire.exon import Exon
from majiq.grimoire.lsv import quant_lsv

from majiq.src.gff import parse_gff3
from majiq.src.constants import *
from voila.splice_graphics import LsvGraphic
import pickle

cdef list accepted_transcripts = ['mRNA', 'transcript']
cdef str transcript_id_keys = 'ID'
cdef list gene_name_keys = ['Name', 'gene_name']
cdef list gene_id_keys = ['ID', 'gene_id']


cdef int _read_gff(str filename, str outDir, object logging=None) except -1:
    """
    :param filename: GFF input filename
    :param list_of_genes: List of genes that will be updated with all the gene_id detected on the gff file
    :param logging: logger object
    :return: :raise RuntimeError:
    """
    cdef list all_genes = []
    cdef dict gene_id_dict = {}
    cdef dict trcpt_id_dict = {}
    cdef dict exon_dict = {}
    cdef dict elem_dict = {}

    cdef str chrom, name, strand
    cdef int start, end
    cdef bint bb
    cdef list ind_list

    for record in parse_gff3(filename):
        chrom = record.seqid
        strand = record.strand
        start = record.start
        end = record.end

        if record.type == 'gene':
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
            all_genes.append((gene_id, gene_name, chrom, strand, start, end))
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
                logging.info("Error, incorrect gff. mRNA %s doesn't have valid gene %s" % (transcript_name, parent))
                raise

        elif record.type == 'exon':
            parent_tx_id = record.attributes['Parent']
            try:
                exon_dict[gn_id].append((start, True))
                exon_dict[gn_id].append((end, False))
                trcpt_id_dict[parent_tx_id][1].append((start, end))

            except KeyError:
                logging.WARNING("Error, incorrect gff. exon at line %s "
                                "doesn't have valid mRNA %s" % (0, parent_tx_id))


    for parent_tx_id, (gn_id, coord_list) in trcpt_id_dict.items():
        last_ss = FIRST_LAST_JUNC
        coord_list.sort(key=lambda x: (x[0], x[1]))
        for xx, yy in coord_list:
            elem_dict[gn_id].append((last_ss, xx, 1, J_TYPE))
            last_ss = yy

        elem_dict[gn_id].append((last_ss, FIRST_LAST_JUNC, 1, J_TYPE))

    merge_exons(exon_dict, elem_dict)
    _dump_elems_list(elem_dict, all_genes, outDir)

    # #test
    # from deleteme import generate_lookuptree
    #
    # generate_lookuptree(all_genes)
    del all_genes


cdef int merge_exons(dict exon_dict, dict elem_dict) except -1:
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
        elem_dict[gne_id] = np.array(elem_dict[gne_id])


#######
# HDF5 API
#######

cdef int _dump_lsv_coverage(out_f, cov_list, attrs_list):
    out_f.create_dataset(JUNCTIONS_DATASET_NAME, data=cov_list,
                         compression='gzip', compression_opts=9)
    out_f.create_dataset(JUNCTIONS_ATTRS, data=attrs_list,
                         compression='gzip', compression_opts=9)


cdef int _dump_elems_list(dict elem_dict, list gene_info, str outDir) except -1:

    dt=np.dtype('|S250, |S250, |S32, S1, u4, u4')
    elem_dict['gene_info'] = np.array(gene_info, dtype=dt)
    np.savez(get_build_temp_db_filename(outDir), **elem_dict)


cdef int _read_junction(np.ndarray row, str gne_id, dict jjs, list exs, list irs, int default_index) except -1:
    jjs[(row[0], row[1])] = Junction(row[0], row[1], gne_id, default_index, annot=bool(row[2]))


cdef int _read_exon(np.ndarray row, str gne_id, dict jjs, list exs, list irs, int default_index) except -1:
    exs.append(Exon(row[0], row[1], annot=bool(row[2])))


cdef int _read_ir(np.ndarray row, str gne_id, dict jjs, list exs, list irs, int default_index) except -1:
    irs.append(Intron(row[0], row[1], annot=bool(row[2]), db_idx=-1))


cdef int _pass_ir(np.ndarray row, str gne_id, dict jjs, list exs, list irs, int default_index) except -1:
    pass


cdef _get_extract_lsv_list(list list_of_lsv_id, list file_list, int msamples):
    cdef list result = []
    cdef int n_exp = len(file_list)
    cdef str lsv_id, lsv_type, fname
    cdef int fidx
    cdef np.ndarray cov, lsv_cov

    for lsv_id in list_of_lsv_id:
        lsv_type = None
        for fidx, fname in enumerate(file_list):

            with h5py.File(fname, 'r') as data:
                try:
                    if lsv_type is None:
                        lsv_type = data['LSVs/%s' % lsv_id].attrs['type']
                        njunc = len(lsv_type.split('|')) -1
                        lsv_cov = np.zeros(shape=(n_exp, njunc, msamples),  dtype=float)

                    assert data['LSVs/%s' % lsv_id].attrs['type'] == lsv_type, "ERROR lsv_type doesn't match " \
                                                                               "for %s" % lsv_id
                    cov = data['LSVs/%s' % lsv_id].attrs['coverage']
                    #lsv_cov[fidx] = data[JUNCTIONS_DATASET_NAME][cov[0]:cov[1]]
                    lsv_cov[fidx] = np.array([data[JUNCTIONS_DATASET_NAME][xidx] for xidx in cov])
                except KeyError:
                    pass

        qq = quant_lsv(lsv_id, lsv_type, lsv_cov)
        result.append(qq)
    return result


####
# API
##

def read_meta_info(list_of_files):
    meta = {'experiments': []}
    m_samples = None
    for fl in list_of_files:
        with h5py.File(fl, 'r') as fp :
            if m_samples is not None:
                assert m_samples == fp.attrs['m_samples'], "uneven number of bootstrap samples"
            else:
                m_samples = fp.attrs['m_samples']
            meta['experiments'].append(fp.attrs['sample_id'])
            try:
                if meta['genome'] != fp.attrs['genome']:
                    raise RuntimeError('Combining experiments from different genome assemblies. Exiting')
            except KeyError:
                meta['genome'] = fp.attrs['genome']
                continue
    meta['m_samples'] = m_samples
    return meta

cpdef extract_lsv_summary(list files, int minnonzero, int min_reads, dict epsi=None, int percent=-1, object logger=None):
    cdef dict lsv_list = {}
    cdef dict lsv_graphic = {}
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
        data = h5py.File(ff, 'r')
        for xx in data['LSVs']:
            mtrx = np.array([data['junc_cov'][xidx] for xidx in data['LSVs/%s' % xx].attrs['coverage']])
            try:
                vals = (mtrx[:, 0] >= min_reads * (mtrx[:, 1] >= minnonzero))
            except IndexError:
                logger.info("Skipping incorrect lsv %s" % xx)
                continue

            try:
                lsv_list[xx][fidx] = int(vals.sum() >= 1)
                if epsi is not None:
                    epsi[xx] += mtrx[:,0]
            except:
                lsv_list[xx]= [0] * nfiles
                lsv_list[xx][fidx] = int(vals.sum() >= 1)
                lsv_graphic[xx] = LsvGraphic.easy_from_hdf5(data['LSVs/%s/visual' % xx])
                if epsi is not None:
                    epsi[xx] = mtrx[:,0]

    if epsi is not None:
        for xx in epsi.keys():
            epsi[xx] = epsi[xx] / nfiles
            epsi[xx] = epsi[xx] / epsi[xx].sum()
            epsi[xx][np.isnan(epsi[xx])] = 1.0 / nfiles

    lsv_id_list = [xx for xx,yy in lsv_list.items() if np.sum(yy) >= percent]

    return lsv_id_list, lsv_graphic


cpdef int dump_lsv_coverage(out_f, cov_list, attrs_list):
    _dump_lsv_coverage(out_f, cov_list, attrs_list)

cpdef int init_majiq_file(str filename, str out_dir, str genome, int msamples):

    with h5py.File('%s/%s.majiq' % (out_dir, filename), 'w') as f:
        f.create_dataset('junctions', (5000, msamples), maxshape=(None, msamples))
        f.create_dataset('junc_cov', (5000, 2), maxshape=(None, 2))

        # fill meta info
        f.attrs['sample_id'] = filename
        f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
        f.attrs['VERSION'] = VERSION
        f.attrs['lsv_idx'] = 0
        f.attrs['num_lsvs'] = 0
        f.attrs['genome'] = genome



cpdef int parse_annot(str filename, str out_dir, object logging=None):

    try:
        _read_gff(filename=filename, outDir=out_dir, logging=logging)
    except Exception as e:
        print e
        raise


cpdef dict retrieve_db_genes(str out_dir):

    cdef list names = ['id', 'name', 'chromosome', 'strand', 'start', 'end']
    cdef dict dd

    gene_info = np.load('%s.npz' % get_build_temp_db_filename(out_dir))['gene_info']
    dd = {gene_info[ii][0].decode('UTF-8'):{xx: gene_info[ii][idx]
                                            for idx, xx in enumerate(names)}
          for ii in range(gene_info.shape[0])}

    return dd

cpdef int retrieve_db(str gne_id, str out_dir, dict dict_junctions,list list_exons, list list_introns,
                     list denovo_ir=[], int default_index=-1):


    cdef dict gne_dict = {}
    cdef dict j_attrs, ex_attrs
    cdef Junction junc
    cdef str xx
    cdef object db_f
    cdef int njuncs = 0


    if list_introns is not None:
        func_list = {EX_TYPE: _read_exon, IR_TYPE: _read_ir, J_TYPE: _read_junction}
    else:
        func_list = {EX_TYPE: _read_exon, IR_TYPE: _pass_ir, J_TYPE: _read_junction}


    all_files = np.load('%s.npz' % get_build_temp_db_filename(out_dir))
    mtrx = all_files[gne_id]
    for i in range(mtrx.shape[0]):
        func_list[mtrx[i,3]](mtrx[i], gne_id, dict_junctions, list_exons, list_introns,
                             default_index)


    return njuncs


cpdef dict retrieve(str gne_id, str out_dir, dict dict_junctions,list list_exons, list list_introns,
                     list denovo_ir=[], int default_index=-1):


    cdef dict gne_dict = {}
    cdef dict j_attrs, ex_attrs
    cdef Junction junc
    cdef str xx
    cdef object db_f
    cdef int njuncs = 0

    cdef list names = ['id', 'name', 'chromosome', 'strand', 'start', 'end']
    cdef dict dd

    gene_info = np.load('%s.npz' % get_build_temp_db_filename(out_dir))['gene_info']




    if list_introns is not None:
        func_list = {EX_TYPE: _read_exon, IR_TYPE: _read_ir, J_TYPE: _read_junction}
    else:
        func_list = {EX_TYPE: _read_exon, IR_TYPE: _pass_ir, J_TYPE: _read_junction}


    all_files = np.load(get_build_temp_db_filename(out_dir))
    gne_list = all_files['gene_info']
    for gne_row in gne_list:
        gne_dict[gne_row[0]] = gne_row
        gne_dict[gne_row[0]] = {xx: gne_row[idx]for idx, xx in enumerate(names)}

        mtrx = all_files[gne_id]
        for i in range(mtrx.shape[0]):
            func_list[mtrx[i,3]](mtrx[i], gne_id, dict_junctions, list_exons, list_introns, default_index)

    return gne_dict



def retrieve_db_info(str gne_id, str out_dir, dict dict_junctions,list list_exons, list list_introns,
                     list denovo_ir=[], int default_index=-1):

    cdef dict j_attrs, ex_attrs
    cdef Junction junc
    cdef str xx
    cdef object db_f
    cdef int njuncs = 0

    if not denovo_ir:
        mode = 'r'
    else:
        mode = 'r+'

    with h5py.File(get_build_temp_db_filename(out_dir), mode=mode) as db_f:

        for xx in db_f['%s/junctions' % gne_id]:
            j_attrs = dict(db_f['%s/junctions/%s' % (gne_id, xx)].attrs)
            njuncs +=1
            dict_junctions[(j_attrs['start'], j_attrs['end'])] = Junction(j_attrs['start'], j_attrs['end'],
                                                                          gne_id, default_index,
                                                                          annot=j_attrs['annotated'])

        mtrx = np.array(db_f['%s/db_coords' % gne_id])

        for idx, row in enumerate(mtrx[mtrx[:, 3] == EX_TYPE, :]):
            list_exons.append(Exon(row[0], row[1], annot=bool(row[2]), db_idx=idx))

        if list_introns is not None:
            for row in mtrx[mtrx[:, 3] == IR_TYPE, :]:
                list_introns.append(Intron(row[0], row[1], annot=bool(row[2]), db_idx=-1))

            if denovo_ir:
                n_mtrx = np.array(denovo_ir)
                shp = mtrx.shape
                shp_new = shp[0] + n_mtrx.shape[0]
                db_f['%s/db_coords' % gne_id].resize((shp_new, shp[1]))
                db_f['%s/db_coords' % gne_id][shp[0]:] = n_mtrx

                for ir in denovo_ir:
                    list_introns.append(Intron(ir[0], ir[1], annot=False, db_idx=-1))


    return njuncs




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

def get_extract_lsv_list(list list_of_lsv_id, list file_list, int msamples):
    return _get_extract_lsv_list(list_of_lsv_id, file_list, msamples)