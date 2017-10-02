import datetime
import gzip
from collections import namedtuple
cimport numpy as np
import h5py
import numpy as np
import urllib.parse as urllib
import os
from majiq.grimoire.junction cimport Junction
from majiq.grimoire.junction import Junction
from majiq.grimoire.exon cimport Exon, Intron
from majiq.grimoire.exon import Exon
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
    cdef dict all_genes = {}
    cdef dict gene_id_dict = {}
    cdef dict trcpt_id_dict = {}
    cdef dict exon_dict = {}

    cdef str chrom, name, strand
    cdef int start, end
    cdef bint bb
    cdef list ind_list
    # cdef set in_juncs = set()

    #with h5py.File(get_build_temp_db_filename(outDir), "w") as db_f:
    db_f = h5py.File(get_build_temp_db_filename(outDir), "w")
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

            _dump_gene(db_f, gene_id, gene_name, chrom, strand, start, end)
            exon_dict[gene_id] = []
            # gn = Gene(gene_id, gene_name, chrom, strand, start, end)
            #exist_overlapping_gene(all_genes[chrom], strand, start, end)

            if gene_id in gene_id_dict:
                raise RuntimeError('Two Different Genes with the same name %s' % gene_name)
            try:
                all_genes[chrom].append([start, end, gene_id])
            except KeyError:
                all_genes[chrom] = [[start, end, gene_id]]
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
                # _dump_junctions(db_f, gn_id, last_ss, start, parent_tx_id, annot=True)
                trcpt_id_dict[parent_tx_id][1].append((start, end))

            except KeyError:
                logging.WARNING("Error, incorrect gff. exon at line %s "
                                "doesn't have valid mRNA %s" % (0, parent_tx_id))

    for parent_tx_id, (gn_id, coord_list) in trcpt_id_dict.items():
        last_ss = FIRST_LAST_JUNC
        coord_list.sort(key=lambda x: (x[0], x[1]))
        for xx, yy in coord_list:
            _dump_junctions(db_f, gn_id, last_ss, xx, parent_tx_id, annot=True)
            last_ss = yy
        _dump_junctions(db_f, gn_id, last_ss, FIRST_LAST_JUNC, parent_tx_id, annot=True)

    merge_exons(db_f, exon_dict)

    db_f.close()
    del all_genes


cdef int merge_exons(db_f, dict exon_dict) except -1:
    cdef list ex_list
    cdef str gne_id
    cdef tuple x
    cdef int ex_start, ex_end, nopen, coord
    cdef bint is_start


    for gne_id, ex_list in exon_dict.items():
        ex_list.sort(key=lambda x:x[0])
        ex_start = -1
        ex_end = -1
        nopen = 0
        for coord, is_start in ex_list:
            #print(coord, is_start, nopen, ex_start, ex_end)
            if is_start:
                if ex_end != -1:
                    #print("NEW EXON, ", gne_id, ex_start, ex_end)
                    _dump_exon(db_f, gne_id, ex_start, ex_end)
                    if nopen > 0:
                        _dump_intron(db_f, gne_id, ex_end+1, coord-1, annot=True)
                    ex_end = -1
                    nopen = 0
                    ex_start = coord

                ex_start = coord if ex_start == -1 or coord < ex_start else ex_start
                nopen += 1

            else:
                nopen -= 1
                ex_end = coord if coord > ex_end else ex_end

        if ex_end != -1:
            #print(gne_id, ex_list)
            _dump_exon(db_f, gne_id, ex_start, ex_end)


#######
# HDF5 API
#######

cdef int _dump_junctions(db_f, str gne_id, int start, int end, str transcript_id, bint annot=False) except -1:
    jid = '%s/junctions/%s-%s' % (gne_id, start, end)

    if jid not in db_f:
        h_jnc = db_f.create_group(jid)
        h_jnc.attrs['start'] = start
        h_jnc.attrs['end'] = end
        h_jnc.attrs['annotated'] = annot
        #h_jnc.attrs['transcript_id_list'] = [transcript_id.encode('utf8')]

cdef int _dump_exon(db_f, str gne_id, int start, int end, bint annot=True) except -1:
    jid = '%s/exons/%s-%s' % (gne_id, start, end)

    if jid not in db_f:
        h_jnc = db_f.create_group(jid)
        h_jnc.attrs['start'] = start
        h_jnc.attrs['end'] = end
        h_jnc.attrs['annot'] = annot


cdef int _dump_gene(db_f, str gne_id, str gne_name, str chrom, str strand, int start, int end) except -1:
    h_gen = db_f.create_group(gne_id)
    h_gen.attrs['id'] = gne_id
    h_gen.attrs['name'] = gne_name
    h_gen.attrs['chromosome'] = chrom
    h_gen.attrs['strand'] = strand
    h_gen.attrs['start'] = start
    h_gen.attrs['end'] = end


cdef int _dump_intron(db_f, str gne_id, int start, int end, bint annot=False) except -1:
    jid = '%s/ir/%s-%s' % (gne_id, start, end)
    if jid not in db_f:
        h_jnc = db_f.create_group(jid)
        h_jnc.attrs['start'] = start
        h_jnc.attrs['end'] = end
        h_jnc.attrs['annotated'] = annot


def junction_to_tmp(gne_id, Junction junc, object hdf5grps):
    cdef str jid = "%s-%s" % (junc.start, junc.end)
    h_jnc = hdf5grps.create_group('%s/%s' % (gne_id, jid))
    h_jnc.attrs['id'] = jid
    h_jnc.attrs['start'] = junc.start
    h_jnc.attrs['end'] = junc.end
    h_jnc.attrs['intronic'] = junc.intronic
    h_jnc.attrs['coverage_index'] = junc.index


####
# API
##

def read_meta_info(list_of_files):
    meta = {'experiments': []}
    for fl in list_of_files:
        with h5py.File(fl, 'r') as fp :
            meta['experiments'].append(fp.attrs['sample_id'])
            try:
                if meta['genome'] != fp.attrs['genome']:
                    raise RuntimeError('Combining experiments from different genome assemblies. Exiting')
            except KeyError:
                meta['genome'] = fp.attrs['genome']
                continue
    return meta

def extract_lsv_summary(list files, int minnonzero, int min_reads, dict epsi=None, int percent=-1, object logger=None):
    cdef dict lsv_list = {}
    cdef dict lsv_graphic = {}
    cdef int nfiles = len(files)
    cdef int fidx
    cdef str ff, xx
    cdef np.ndarray mtrx, vals

    if percent == -1:
        percent = nfiles / 2
        percent = percent + 1 if nfiles % 2 != 0 else percent


    for fidx, ff in enumerate(files):
        if logger:
            logger.info("Parsing file: %s" % ff)
        data = h5py.File(ff, 'r')
        for xx in data['LSVs']:
            mtrx = data['junc_cov'][data['LSVs/%s' % xx].attrs['coverage'][0]:data['LSVs/%s' % xx].attrs['coverage'][1]]
            vals = (mtrx[:, 0] >= min_reads * (mtrx[:, 1] >= minnonzero))
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
            epsi[xx] /= nfiles
            epsi[xx] /= (epsi[xx].sum() - epsi[xx])
            epsi[xx][np.isnan(epsi[xx])] = 0.5

    lsv_id_list = [xx for xx,yy in lsv_list.items() if np.sum(yy) > percent]

    return lsv_id_list, lsv_graphic


def dump_junctions(db_f, str gne_id, int start, int end, str transcript_id, bint annot=False):
    _dump_junctions(db_f, gne_id, start, end, transcript_id, annot=annot)

def dump_intron(db_f, str gne_id, int start, int end, bint annot=False):
    _dump_intron(db_f, gne_id, start, end, annot=annot)

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

cpdef retrieve_db_genes(str out_dir):

    cdef dict dd = {}
    with h5py.File(get_build_temp_db_filename(out_dir), 'r') as db_f:
        for ii in db_f:
            dd[ii] = dict(db_f[ii].attrs)
    return dd


cpdef list get_list_of_genes(str out_dir):
    cdef list list_of_genes
    with h5py.File(get_build_temp_db_filename(out_dir), 'r') as db_f:
        list_of_genes = list(db_f.keys())
    return list_of_genes


cpdef int parse_annot(str filename, str out_dir, object logging=None):

    try:
        _read_gff(filename=filename, outDir=out_dir, logging=logging)
    except Exception as e:
        print e
        raise


cpdef np.ndarray get_covered_junctions(str gne_id, dict dict_junctions, list list_exons, list fitfunc_r, list sam_list,
                                       int readLen, str out_dir):

    cdef int nexperiments = len(sam_list)
    cdef dict dd_jj = {}
    cdef int exp_idx, njunc, idx
    cdef str jjid
    cdef Junction junc
    cdef dict junc_attrs
    cdef np.ndarray junc_mtrx
    cdef int eff_readsize = (readLen - 2*MIN_BP_OVERLAP) + 1

    for exp_idx in range(nexperiments):
        with h5py.File(get_builder_temp_majiq_filename(out_dir, sam_list[exp_idx]), 'r') as fp:
            fitfunc_r.append(fp.attrs['one_over_r'])
            try:
                for jjid in fp[gne_id]:
                    junc_attrs = dict(fp['%s/%s' % (gne_id, jjid)].attrs)
                    try:
                        dd_jj[jjid][exp_idx] = junc_attrs['coverage_index']
                    except KeyError:
                        junc = Junction(junc_attrs['start'], junc_attrs['end'], gne_id, -1)
                        dict_junctions[(junc_attrs['start'], junc_attrs['end'])] = junc
                        dd_jj[jjid] = [-1] * nexperiments
                        dd_jj[jjid][exp_idx] = junc_attrs['coverage_index']
            except KeyError as e:
                # print('ERROR', e, jjid, junc_attrs)
                continue

    njunc = len(dd_jj)
    junc_mtrx = np.zeros(shape=(nexperiments, njunc, eff_readsize), dtype=np.float)

    for exp_idx in range(nexperiments):
        with h5py.File(get_builder_temp_majiq_filename(out_dir, sam_list[exp_idx]), 'r') as fp:
            for idx, jjid in enumerate(dd_jj.keys()):
                if dd_jj[jjid][exp_idx] > -1:
                    junc_mtrx[exp_idx, idx,] = fp[JUNCTIONS_DATASET_NAME][dd_jj[jjid][exp_idx]]
                dict_junctions[tuple([int(xx) for xx in jjid.split('-')])].index = idx

    del dd_jj
    return junc_mtrx


def retrieve_db_info(str gne_id, str out_dir, list list_exons, dict dict_junctions, list list_introns):

    cdef dict j_attrs, ex_attrs
    cdef Junction junc
    cdef str xx
    cdef object db_f
    cdef int njuncs = 0

    with h5py.File(get_build_temp_db_filename(out_dir), 'r') as db_f:
        if dict_junctions is not None:
            for xx in db_f['%s/junctions' % gne_id]:
                j_attrs = dict(db_f['%s/junctions/%s' % (gne_id, xx)].attrs)
                njuncs +=1
                dict_junctions[(j_attrs['start'], j_attrs['end'])] = Junction(j_attrs['start'], j_attrs['end'],
                                                                              gne_id, -1, annot=j_attrs['annotated'])

        for xx in db_f['%s/exons' % gne_id]:
            ex_attrs = dict(db_f['%s/exons/%s' % (gne_id, xx)].attrs)
            list_exons.append(Exon(ex_attrs['start'], ex_attrs['end'], annot=ex_attrs['annot']))
        list_exons.sort(key=lambda xx: (xx.start, xx.end))

        if list_introns is not None and 'ir' in db_f['%s' % gne_id]:
            for xx in db_f['%s/ir' % gne_id]:
                ir_attrs = dict(db_f['%s/ir/%s' % (gne_id, xx)].attrs)
                list_introns.append(Intron(ir_attrs['start'], ir_attrs['end'], annot=ir_attrs['annotated']))
            list_introns.sort(key=lambda xx: (xx.start, xx.end))

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