import datetime
import gzip
from collections import namedtuple
cimport numpy as np
import h5py
import numpy as np
import urllib.parse as urllib
from majiq.grimoire.junction cimport Junction
from majiq.grimoire.junction import Junction
from majiq.grimoire.exon cimport Exon
from majiq.grimoire.exon import Exon
from majiq.src.gff import parse_gff3
from majiq.src.constants import *
from majiq.src.io_bam cimport read_juncs_from_bam

# import majiq.src.utils as majiq_utils

#
# cdef list gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
# GFFRecord = namedtuple("GFFRecord", gffInfoFields)
#
#
# cdef dict __parse_gff_attributes(str attribute_string):
#     """
#     Parse the GFF3 attribute column and return a dict
#     :param attribute_string:
#     """
#     cdef dict ret
#     cdef str key, value
#
#     if attribute_string == ".":
#         return {}
#     ret = {}
#     for attribute in attribute_string.split(";"):
#         key, value = attribute.split("=")
#         key = urllib.unquote(key)
#         if key in ret:
#             key = 'extra_%s' % key
#             if key not in ret:
#                 ret[key] = []
#             ret[key].append(urllib.unquote(value))
#         else:
#             ret[key] = urllib.unquote(value)
#     return ret


# def __parse_gff3(str filename):
#     """
#     A minimalistic GFF3 format parser.
#     Yields objects that contain info about a single GFF3 feature.
#
#     Supports transparent gzip decompression.
#     """
#     # Parse with transparent decompression
#
#     # cdef object infile
#     cdef str line
#     # cdef list parts
#     cdef dict normalized_info
#     open_func = gzip.open if filename.endswith(".gz") else open
#     with open_func(filename) as infile:
#
#         for line in infile:
#
#             if line.startswith("#"):
#                 continue
#             parts = line.strip().split("\t")
#             # If this fails, the file format is not standard-compatible
#             assert len(parts) == len(gffInfoFields)
#             # Normalize data
#             normalized_info = {
#                 "seqid": None if parts[0] == "." else urllib.unquote(parts[0]),
#                 "source": None if parts[1] == "." else urllib.unquote(parts[0]),
#                 "type": None if parts[2] == "." else urllib.unquote(parts[2]),
#                 "start": None if parts[3] == "." else int(parts[3]),
#                 "end": None if parts[4] == "." else int(parts[4]),
#                 "score": None if parts[5] == "." else float(parts[5]),
#                 "strand": None if parts[6] == "." else urllib.unquote(parts[6]),
#                 "phase": None if parts[7] == "." else urllib.unquote(parts[7]),
#                 "attributes": __parse_gff_attributes(parts[8])
#             }
#             # Alternatively, you can emit the dictionary here, if you need mutabwility:
#             yield normalized_info
#             # yield GFFRecord(**normalized_info)




cdef inline int __read_STAR_junc_file(str filename, set in_juncs, list out_jj, dict out_dd) except -1:

    with open(filename, 'r') as fp:
        for ln in fp.readlines():
            tab = ln.strip().split()
            if (tab[0], tab[1], tab[2]) in in_juncs:
                continue
            try:
                out_dd[(tab[0], tab[1], tab[2])] += tab[6]
            except KeyError:
                out_dd[(tab[0], tab[1], tab[2])] = tab[6]
            out_jj.append((tab[0], tab[1], tab[2], tab[3]))


cdef int _read_juncs(db_f, list file_list, set in_juncs, dict all_genes, int min_denovo) except -1:

    cdef str ln, fname, chrom='', jjid
    cdef list tab
    cdef list jj_list = []
    cdef dict jj_dd = {}
    cdef tuple jj
    cdef int gidx=0, ngenes=0
    cdef file fp

    for bb, fname in file_list:
        if bb:
            __read_STAR_junc_file(fname, in_juncs, jj_list, jj_dd)
        else:
            read_juncs_from_bam(fname, in_juncs, jj_list, jj_dd)

    jj_list.sort(key=lambda xx: (xx[0], xx[1], xx[2]))
    jj_dd = {kk: float(vv)/len(file_list) for kk, vv in jj_dd.items()}

    for jj in jj_list:
        if jj_dd[(jj[0], jj[1], jj[2])] < min_denovo:
            continue

        if chrom != jj[0]:
            chrom = jj[0]
            ngenes = len(all_genes[chrom])

        while gidx < ngenes :
            gstart, gend, gid = all_genes[chrom][gidx]
            if gend < jj[1]:
                gidx += 1
                continue
            if gstart > jj[2]:
                break

            if gstart<=jj[1] and gend>= jj[2]:
                jjid = '%s/junctions/%s-%s' % (gid, jj[1], jj[2])
                if jjid not in db_f:
                    _dump_junctions(db_f, gid, jj[1], jj[2], None, annot=False)
                break


cdef list accepted_transcripts = ['mRNA', 'transcript']
cdef str transcript_id_keys = 'ID'
cdef list gene_name_keys = ['Name', 'gene_name']
cdef list gene_id_keys = ['ID', 'gene_id']


cdef int _read_gff(str filename, str outDir, dict exp_groups, bint denovo, int min_denovo,
                   object logging=None) except -1:
    """
    :param filename: GFF input filename
    :param list_of_genes: List of genes that will be updated with all the gene_id detected on the gff file
    :param logging: logger object
    :return: :raise RuntimeError:
    """
    cdef dict all_genes = {}
    cdef dict gene_id_dict = {}
    cdef dict trcpt_id_dict = {}
    cdef list exon_list = []

    cdef str chrom, name, strand
    cdef int start, end
    cdef bint bb
    cdef list ind_list

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
                trcpt_id_dict[record.attributes['ID']] = [gn_id, 0]

            except KeyError:
                logging.info("Error, incorrect gff. mRNA %s doesn't have valid gene %s" % (transcript_name, parent))
                raise

        elif record.type == 'exon':
            parent_tx_id = record.attributes['Parent']
            try:
                gn_id, last_ss = trcpt_id_dict[parent_tx_id]
                ex_tuple = (start, (start, end), False)
                exon_list.append(ex_tuple)
                ex_tuple = (end, (start, end), True)
                exon_list.append(ex_tuple)
                _dump_exon(db_f, gn_id, start, end)

                if last_ss > 0:
                    _dump_junctions(db_f, gn_id, last_ss, start, parent_tx_id, annot=True)
                trcpt_id_dict[parent_tx_id][1] = end

            except KeyError:
                logging.WARNING("Error, incorrect gff. exon at line %s "
                                "doesn't have valid mRNA %s" % (0, parent_tx_id))

    if denovo:

        for name, ind_list in exp_groups.items():
            _read_juncs(db_f, ind_list, all_genes, min_denovo)

    db_f.close()
    del all_genes

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
    else:
        h_jnc = db_f[jid]
        #h_jnc.attrs['transcript_id_list'].append(transcript_id.encode('utf8'))

cdef int _dump_exon(db_f, str gne_id, int start, int end) except -1:
    jid = '%s/exons/%s-%s' % (gne_id, start, end)

    if jid not in db_f:
        h_jnc = db_f.create_group(jid)
        h_jnc.attrs['start'] = start
        h_jnc.attrs['end'] = end


cdef int _dump_gene(db_f, str gne_id, str gne_name, str chrom, str strand, int start, int end) except -1:
    h_gen = db_f.create_group(gne_id)
    h_gen.attrs['id'] = gne_id
    h_gen.attrs['name'] = gne_name
    h_gen.attrs['chromosome'] = chrom
    h_gen.attrs['strand'] = strand
    h_gen.attrs['start'] = start
    h_gen.attrs['end'] = end


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


cpdef int parse_annot(str filename, str out_dir, dict exp_groups, bint denovo, int min_denovo, object logging=None):

    _read_gff(filename=filename, outDir=out_dir, exp_groups=exp_groups, denovo=denovo, min_denovo=min_denovo,
              logging=logging)


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


cpdef int retrieve_db_info(str gne_id, str out_dir, list list_exons, dict dict_junctions):

    cdef dict j_attrs, ex_attrs
    cdef Junction junc
    cdef str xx
    cdef object db_f

    with h5py.File(get_build_temp_db_filename(out_dir), 'r') as db_f:
        if dict_junctions is not None:
            for xx in db_f['%s/junctions' % gne_id]:
                j_attrs = dict(db_f['%s/junctions/%s' % (gne_id, xx)].attrs)
                dict_junctions[(j_attrs['start'], j_attrs['end'])] = Junction(j_attrs['start'], j_attrs['end'],
                                                                              gne_id, -1, annot=True)

        for xx in db_f['%s/exons' % gne_id]:
            ex_attrs = dict(db_f['%s/exons/%s' % (gne_id, xx)].attrs)
            list_exons.append(Exon(ex_attrs['start'], ex_attrs['end'], annot=True))
        list_exons.sort(key=lambda xx: (xx.start, xx.end))