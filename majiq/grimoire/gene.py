#!/usr/bin/python
import gc

from majiq.grimoire.exon import Exon, ExonTx, collapse_list_exons, detect_exons
from majiq.grimoire.junction import Junction
from majiq.grimoire.lsv import LSV
#from majiq.src import config_old as majiq_config
from majiq.src.config import Config
from majiq.grimoire.exon import new_exon_definition
from majiq.src.constants import *


class Gene:

    def __init__(self, gene_id, gene_name, chrom, strand, start, end, retrieve=False):
        self.id = gene_id
        self.name = gene_name
        self.chromosome = chrom
        self.strand = strand
        self.exons = []
        self.start = start
        self.end = end
        self.ir_definition = []
        self.antis_gene = []

        self.total_read = 0
        if not retrieve:
            self.transcript_tlb = {}
            self.temp_txex_list = []
        self.junc_index = 0
        self.gc_content = None

    def __del__(self):
        if hasattr(self, "junc_cov"):
            del self.junc_cov
        if hasattr(self, "junc_pos"):
            del self.junc_pos
        if hasattr(self, "gc_content"):
            del self.gc_content

        for ex in self.exons:
            del ex

    def __hash__(self):
        return hash((self.id, self.chromosome, self.strand, self.start, self.end))

    def to_hdf5(self, hdf5grps):
        h_gen = hdf5grps.create_group("%s" % self.id)

        h_gen.attrs['id'] = self.id
        h_gen.attrs['name'] = self.name
        h_gen.attrs['chromosome'] = self.chromosome
        h_gen.attrs['strand'] = self.strand
        h_gen.attrs['start'] = self.start
        h_gen.attrs['end'] = self.end
        if len(self.ir_definition) > 0:
            h_gen.attrs['ir_definition'] = self.ir_definition
        if len(self.antis_gene) > 0:
            h_gen.attrs['antis_gene'] = self.antis_gene

        [ex.to_hdf5(h_gen) for ex_count, ex in enumerate(self.exons)]

    @staticmethod
    def get_junctions_from_hdf5(hdf5_gene):
        junc_res = []
        for ex_hdf5 in hdf5_gene['exons']:
            for extx_hdf5 in hdf5_gene['exons/%s' % ex_hdf5]:
                for jj in hdf5_gene['exons/%s/%s/p5_junc' % (ex_hdf5, extx_hdf5)]:
                    junc = hdf5_gene['exons/%s/%s/p5_junc/%s' % (ex_hdf5, extx_hdf5, jj)]
                    junc_res.append((junc.attrs['start'], junc.attrs['end']))

        return junc_res

    @staticmethod
    def get_exons_from_hdf5(hdf5_gene):
        ex_res = []
        for ex_hdf5 in hdf5_gene['exons']:
            exn = hdf5_gene['exons/%s' % ex_hdf5]
            ex_res.append((exn.attrs['start'], exn.attrs['end']))

        return ex_res

    def get_id(self):
        return self.id

    def get_name(self):
        return self.name

    def get_strand(self):
        return self.strand

    def get_read_count(self):
        return self.total_read

    def get_exon_list(self):
        return self.exons

    def get_chromosome(self):
        return self.chromosome

    def get_coordinates(self):
        return self.start, self.end

    def get_transcript(self, trans_id):
        return self.transcript_tlb[trans_id]

    def get_ir_definition(self):
        return self.ir_definition

    def incr_junc_index(self):
        self.junc_index += 1

    def get_junc_index(self):
        return self.junc_index

    def set_junc_index(self, val):
        self.junc_index = val

    def get_all_junctions(self, filter=True):
        lst = set()
        for ex in self.get_exon_list():
            for ex_rna in ex.exonRead_list:
                if len(ex_rna.p3_junc) > 0:
                    lst = lst.union(set(ex_rna.p3_junc))
                if len(ex_rna.p5_junc) > 0:
                    lst = lst.union(set(ex_rna.p5_junc))
            for ex_tx in ex.exonTx_list:
                if len(ex_tx.p3_junc) > 0:
                    lst = lst.union(set(ex_tx.p3_junc))
                if len(ex_tx.p5_junc) > 0:
                    lst = lst.union(set(ex_tx.p5_junc))

        s_junc = list(lst)
        return sorted([xx for xx in s_junc if not xx is None and (not filter or xx.get_coverage().sum() > 0)])

    def get_all_introns(self):
        lintrons = []
        if len(self.exons) > 1:
            for ex_idx, ex in enumerate(self.exons[:-1]):
                lintrons.append((ex, self.exons[ex_idx+1]))
        return lintrons

    def add_exon(self, exon):
        self.exons.append(exon)

    def add_transcript(self, tcrpt):
        if tcrpt.txstart < self.start:
            self.start = tcrpt.txstart
        if tcrpt.txend > self.end:
            self.end = tcrpt.txend

        self.transcript_tlb[tcrpt.get_id()] = tcrpt

    def add_ir_definition(self, start, end):
        self.ir_definition.append((start, end))

    def add_read_count(self, read_num):
        self.total_read += read_num

    def exist_antisense_gene(self, list_of_genes):
        # strnd = '-'
        # if self.strand == '-':
        #     strnd = '+'
        for strnd in ('+', '-'):
            for gg in list_of_genes[strnd]:
                if gg.get_id() == self.id:
                    continue
                if self.overlaps(gg):
                    # coords = gg.get_coordinates()
                    # if self.start < coords[1] and self.end > coords[0]:
                    self.set_antisense_gene(gg.get_id())
                    gg.set_antisense_gene(self.id)

    def set_antisense_gene(self, gn_id):
        self.antis_gene.append(gn_id.encode('utf8'))

    def overlaps(self, gne):
        return self.start < gne.end and self.end > gne.start

    def new_annotated_exon(self, start, end, transcript, bl=True, intron=False):
        for txex in self.temp_txex_list:
            if txex.start == start and txex.end == end:
                res = txex
                break
        else:
            res = ExonTx(start, end, transcript, intron=intron)
            if bl:
                self.temp_txex_list.append(res)
        return res

    def exist_junction(self, start, end):
        if start is None or end is None:
            return
        res = None
        for txcpt in self.transcript_tlb.values():
            res = txcpt.in_junction_list(start, end)
            if res is not None:
                break
        return res

    def new_annotated_junctions(self, start, end, trcpt):
        junc = self.exist_junction(start, end)
        if junc is None:
            junc = Junction(start, end, None, None, self.get_id(), annotated=True)
        junc.add_transcript(trcpt)

        return junc

    def prepare_exons(self):
        self.exons.sort()

    def remove_temp_attributes(self):
        del self.temp_txex_list

    def collapse_exons(self):

        self.temp_txex_list.sort()
        list_ex = collapse_list_exons(self.temp_txex_list, self)
        self.exons.extend(list_ex)
        self.prepare_exons()
        self.remove_temp_attributes()

    def get_all_ss(self):
        ss = set()
        for ex in self.exons:
            tx_list = ex.get_annotated_exon()
            for txex in tx_list:
                for junc in txex.get_3prime_junc():
                    coord = junc.end
                    if coord is not None:
                        ss.add((coord, '3prime', junc))
                for junc in txex.get_5prime_junc():
                    coord = junc.start
                    if coord is not None:
                        ss.add((coord, '5prime', junc))

        return sorted(ss)

    def get_exon_by_id(self, ex_id):
        # return self.exons[ex_id-1]
        for ex in self.exons:
            if ex.get_id() == ex_id:
                res = ex
                break
        else:
            res = None

        return res

    def exist_exon(self, start, end):
        """
         .. function: exist_exon (self, start, end):
            Check if the pair (start, end) are in a known exon in the gene. If not return None. We assume
            that the given exon is for the same chromosome and strand than the gene.

            :param start: Start position of the input exon.
            :param end: End Position of te input exon.
            :rtype: exon instance or None
        """

        res = None

        for ee in self.exons:
            ex_str, ex_end = ee.get_coordinates()
            if start < ex_end and end > ex_str:
                res = ee
                break

        return res

    def check_antisense_junctions_hdf5(self, jstart, jend, h5_file):
        majiq_config = Config()
        res = False
        for anti_g in self.antis_gene:
            gg = h5_file[anti_g]
            gg_start = gg.attrs['start']
            gg_end = gg.attrs['end']
            if jend < gg_start or jstart > gg_end:
                continue
            j_list = Gene.get_junctions_from_hdf5(gg)
            for j_st, j_ed in j_list:
                if j_st > jstart or (j_st == jstart and j_ed > jend):
                    break
                elif jstart == j_st and jend == j_ed:
                    res = True
                    break
            if not res:
                for ex_start, ex_end in Gene.get_exons_from_hdf5(gg):

                    coords = [ex_start - majiq_config.get_max_denovo_difference(),
                              ex_end + majiq_config.get_max_denovo_difference()]

                    if coords[0] <= jend <= coords[1] or coords[0] <= jstart <= coords[1]:
                        res = True
                        break
            else:
                break
        return res

    def new_lsv_definition(self, exon, jlist, lsv_type):

        coords = exon.get_coordinates()
        lsv_id = "%s:%d-%d:%s" % (self.get_id(), coords[0], coords[1], lsv_type)
        return LSV(exon, lsv_id, jlist, lsv_type)

    def reset_to_db(self):
        for xx in self.get_all_junctions():
            xx.coverage.fill(0)

    def simplify(self):
        majiq_config = Config()
        jj_set = set()
        for ex in self.exons:
            jlist = []
            for ttype in ('5prime', '3prime'):
                for xx in ex.get_junctions(ttype):
                    if xx is None:
                        continue
                    if xx.get_donor() is None or xx.get_acceptor() is None:
                        jj_set.add(xx)
                    else:
                        jlist.append(xx)

                for exp_idx in range(majiq_config.num_experiments):
                    cover = [float(junc.get_coverage_sum(exp_idx)) for junc in jlist]
                    if sum(cover) == 0:
                        continue
                    bool_map = [majiq_config.simplify_type == SIMPLIFY_ALL or
                                (junc.is_annotated() and majiq_config.simplify_type == SIMPLIFY_DB) or
                                (not junc.is_annotated() and majiq_config.simplify_type == SIMPLIFY_DENOVO)
                                for junc in jlist]

                    jj_set = jj_set.union(set([junc for eidx, junc in enumerate(jlist)
                                               if cover[eidx]/sum(cover) >= majiq_config.simplify_threshold and
                                               bool_map[eidx]]))

            del ex
        self.exons = []
        splice_list = set()
        for jj in jj_set:
            jj.add_donor(None)
            jj.add_acceptor(None)
            splice_list.add((jj.start, '5prime', jj))
            splice_list.add((jj.end, '3prime', jj))
        detect_exons(self, list(splice_list), retrieve=True)
        del splice_list


class Transcript(object):

    def __init__(self, name, gene, txstart, txend):
        self.id = name
        self.gene_id = gene.get_id()
        self.exon_list = []
        self.txstart = txstart
        self.txend = txend

    def get_gene(self):
        majiq_config = Config()
        return majiq_config.gene_tlb[self.gene_id]

    def get_exon_list(self):
        return self.exon_list

    def get_id(self):
        return self.id

    def get_junction_list(self):
        res = set()
        for ex in self.exon_list:
            res.update(set(ex.get_5prime_junc()))
            res.update(set(ex.get_3prime_junc()))

        return sorted(res)

    def prepare_exon_list(self):
        self.exon_list.sort()
        return self.exon_list

    def in_junction_list(self, start, end):
        res = None
        for ex in self.exon_list:
            for jj in ex.get_5prime_junc():
                jjstart, jjend = jj.get_coordinates()
                if jjstart == start and jjend == end:
                    res = jj
                    break
            else:
                continue
            break
        return res

    def add_exon(self, exon):

        self.exon_list.append(exon)
        return


def recreate_gene_tlb(gene_list):
    majiq_config = Config()
    for gn in gene_list:
        majiq_config.gene_tlb[gn.get_id()] = gn


def clear_gene_tlb():
    majiq_config = Config()
    majiq_config.gene_tlb.clear()
    gc.collect()


def retrieve_gene(gene_id, dbfile, all_exp=False, junction_list=None, logger=None):
    majiq_config = Config()
    gg = dbfile[gene_id]
    gg_attrs = dict(gg.attrs)
    try:
        gn = Gene(gene_id, gg_attrs['name'], gg_attrs['chromosome'], gg_attrs['strand'],
                  gg_attrs['start'], gg_attrs['end'], retrieve=True)

        majiq_config.gene_tlb[gene_id] = gn

        if junction_list is None:
            junction_list = {}

        num_exp = 0 if not all_exp else majiq_config.num_experiments

        for ex_grp_id in gg['exons']:
            ex_grp = gg['exons/%s' % ex_grp_id]
            ex_grp_attrs = dict(ex_grp.attrs)
            ex = Exon(ex_grp_attrs['start'], ex_grp_attrs['end'], gn.get_id(),
                      annot=True, isintron=False, indata=ex_grp_attrs['in_data'], retrieve=True)
            gn.exons.append(ex)
            try:
                ex.set_pcr_score(ex_grp_attrs['pcr_name'], ex_grp_attrs['score'], ex_grp_attrs['candidate'])
            except KeyError:
                pass

            if majiq_config.gcnorm:
                ex.set_gc_content_val(ex_grp_attrs['gc_content'])

            for ex_tx_id in ex_grp['tx']:
                ex_tx = ex_grp['tx/%s' % ex_tx_id]
                ex_tx_attrs = dict(ex_tx.attrs)
    #            TODO: Do we need trasncript? for now is None

                transcript_id = None
                ext = ExonTx(ex_tx_attrs['start'], ex_tx_attrs['end'], transcript_id, intron=False)
                ext.gene_name = gene_id
                ex.add_exon_tx(ext)
                ex.add_ss_3p(ex_tx_attrs['start'])
                ex.add_ss_5p(ex_tx_attrs['end'])
                for jj_grp_id in ex_tx["p3_junc"]:
                    jj_grp_attrs = dict(ex_tx["p3_junc/%s" % jj_grp_id].attrs)
                    try:
                        junc = junction_list[jj_grp_attrs['start'], jj_grp_attrs['end']]
                    except KeyError:
                        junc = Junction(jj_grp_attrs['start'], jj_grp_attrs['end'], None, None,
                                        gene_id, annotated=True, retrieve=True, num_exp=num_exp, jindex=-1)
                        junc.donor_id = jj_grp_attrs['donor_id']
                        junc.acceptor_id = jj_grp_attrs['acceptor_id']
                        junc.transcript_id_list = jj_grp_attrs['transcript_id_list']
                        junction_list[jj_grp_attrs['start'], jj_grp_attrs['end']] = junc

                    ext.add_3prime_junc(junc)

                for jj_grp_id in ex_tx["p5_junc"]:
                    jj_grp_attrs = dict(ex_tx["p5_junc/%s" % jj_grp_id].attrs)
                    try:
                        junc = junction_list[jj_grp_attrs['start'], jj_grp_attrs['end']]
                    except KeyError:
                        junc = Junction(jj_grp_attrs['start'], jj_grp_attrs['end'], None, None,
                                        gene_id, annotated=True, retrieve=True, num_exp=num_exp)
                        junc.donor_id = jj_grp_attrs['donor_id']
                        junc.acceptor_id = jj_grp_attrs['acceptor_id']
                        junction_list[jj_grp_attrs['start'], jj_grp_attrs['end']] = junc

                    ext.add_5prime_junc(junc)
    except KeyError:
        logger.info('ERROR in Gene %s: Annotation db analysis is corrupted' % gene_id)
        raise
    return gn


def extract_junctions_hdf5(gene_obj, jj_grp, junction_list, annotated=True, all_exp=False):
    majiq_config = Config()
    num_exp = 1 if not all_exp else majiq_config.num_experiments

    try:
        junc = junction_list[(jj_grp.attrs['start'], jj_grp.attrs['end'])]
        if junc.get_index() == -1:
            junc.idx = gene_obj.get_junc_index()
            gene_obj.incr_junc_index()

    except KeyError:
        if jj_grp.attrs['end'] - jj_grp.attrs['start'] == 1:
            intronic = True
        else:
            intronic = False
        junc = Junction(jj_grp.attrs['start'], jj_grp.attrs['end'], None, None,
                        gene_obj.get_id(), annotated=annotated, retrieve=True, num_exp=num_exp,
                        jindex=gene_obj.get_junc_index(), intronic=intronic)
        junction_list[(jj_grp.attrs['start'], jj_grp.attrs['end'])] = junc
        gene_obj.incr_junc_index()

    return junc


def find_intron_retention(gene_obj, dict_of_junctions, nondenovo, logging=None):
    intron_list = gene_obj.get_all_introns()
    for exon1, exon2 in intron_list:
        ex1_end = exon1.get_coordinates()[1]
        ex2_start = exon2.get_coordinates()[0]
        intron_start = ex1_end + 1
        intron_end = ex2_start - 1

        intron_len = intron_end - intron_start
        if intron_len <= 0:
            continue

        try:
            jin = dict_of_junctions[intron_start]
        except KeyError:
            jin = None
        try:
            jout = dict_of_junctions[intron_end]
        except KeyError:
            jout = None

        if jin is None or jout is None:
            continue
        else:

            exnum = new_exon_definition(intron_start, intron_end, jin, jout, gene_obj, nondenovo=nondenovo,
                                        isintron=True)
            if exnum == -1:
                continue
            logging.debug("NEW INTRON RETENTION EVENT %s, %d-%d" % (gene_obj.get_name(), intron_start, intron_end))
            jin.add_donor(exon1)
            for ex in exon1.exonRead_list:
                st, end = ex.get_coordinates()
                if end == jin.get_coordinates()[0]:
                    ex.add_5prime_junc(jin)
                    break

            jout.add_acceptor(exon2)
            for ex in exon2.exonRead_list:
                st, end = ex.get_coordinates()
                if st == jout.get_coordinates()[1]:
                    ex.add_3prime_junc(jout)
                    break
