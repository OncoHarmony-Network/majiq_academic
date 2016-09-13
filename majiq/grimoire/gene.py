#!/usr/bin/python
import gc

from majiq.grimoire.exon import Exon

from exon import ExonTx, collapse_list_exons
from majiq.grimoire.junction import Junction
from majiq.grimoire.lsv import LSV
from majiq.src import config as majiq_config


class Gene:
    __eq__ = lambda self, other: (self.chromosome == other.chromosome and self.strand == other.strand
                                  and self.start < other.end and self.end > other.start
                                  and self.strand == other.strand)
    __ne__ = lambda self, other: (self.chromosome != other.chromosome or self.strand != other.strand
                                  or self.start >= other.end or self.end <= other.start)
    __lt__ = lambda self, other: (self.chromosome < other.chromosome
                                  or (self.chromosome == other.chromosome
                                      and (self.end < other.start
                                           or (self.end > other.start and self.start < other.end
                                               and self.strand == '+' and other.strand == '-'))))
    __gt__ = lambda self, other: (self.chromosome > other.chromosome
                                  or (self.chromosome == other.chromosome
                                      and (self.start > other.end
                                           or (self.end > other.start and self.start < other.END
                                               and self.strand == '-' and other.strand == '+'))))

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

        #self.readNum = np.zeros(shape=majiq_config.num_experiments, dtype=np.int)
        self.total_read = 0
        if not retrieve:
            self.transcript_tlb = {}
            self.temp_txex_list = []
        else:
            self.lsv_list = []

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
            for extx_hdf5 in hdf5_gene['exons/%s' % ex_hdf5]:
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

    def get_all_junctions(self):
        lst = set()
        for ex in self.get_exon_list():
            for ex_rna in ex.exonRead_list:
                if len(ex_rna.p5_junc) > 0:
                    lst = lst.union(set(ex_rna.p5_junc))
            for ex_tx in ex.exonTx_list:
                if len(ex_tx.p5_junc) > 0:
                    lst = lst.union(set(ex_tx.p5_junc))

        s_junc = list(lst)
        return sorted([xx for xx in s_junc if not xx is None])

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
        self.antis_gene.append(gn_id)

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
            junc = Junction(start, end, None, None, self, annotated=True)
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
                    coord = junc.get_ss_3p()
                    if coord is not None:
                        ss.add((coord, '3prime', junc))
                for junc in txex.get_5prime_junc():
                    coord = junc.get_ss_5p()
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


            # retrieve_gene(anti_g, majiq_config.dbfile)
            # # if not self.antis_gene is None:
            # gg = majiq_config.gene_tlb[anti_g]
            #
            # cc = gg.get_coordinates()
            # if jend < cc[0] or jstart > cc[1]:
            #     continue
            # j_list = gg.get_all_junctions()
            # for jj in j_list:
            #     if not jj.is_annotated():
            #         continue
            #     (j_st, j_ed) = jj.get_coordinates()
            #     if j_st > jstart or (j_st == jstart and j_ed > jend):
            #         break
            #     elif jstart == j_st and jend == j_ed:
            #         res = True
            #         break
            # if not res:
            #     for ex in gg.get_exon_list():
            #         coords = ex.get_coordinates()
            #         # if jend > coords[0]:
            #         #     break
            #         coords = [coords[0] - majiq_config.get_max_denovo_difference(),
            #                   coords[1] + majiq_config.get_max_denovo_difference()]
            #
            #         if coords[0] <= jend <= coords[1] or coords[0] <= jstart <= coords[1]:
            #             res = True
            #             break
            #     del majiq_config.gene_tlb[anti_g]
            # else:
            #     del majiq_config.gene_tlb[anti_g]
            #     break

        return res

    def new_lsv_definition(self, exon, jlist, lsv_type):

        coords = exon.get_coordinates()
        lsv_id = "%s:%d-%d:%s" % (self.get_id(), coords[0], coords[1], lsv_type)
        for lsv in self.lsv_list:
            if lsv.id == lsv_id:
                ret = lsv
                break
        else:
            ret = LSV(exon, lsv_id, jlist, lsv_type)
            self.lsv_list.append(ret)
        return ret


class Transcript(object):

    def __init__(self, name, gene, txstart, txend):
        self.id = name
        self.gene_id = gene.get_id()
        self.exon_list = []
        self.txstart = txstart
        self.txend = txend

    def get_gene(self):
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
    for gn in gene_list:
        majiq_config.gene_tlb[gn.get_id()] = gn


def clear_gene_tlb():
    majiq_config.gene_tlb.clear()
    gc.collect()


def retrieve_gene(gene_id, dbfile, all_exp=False, logger=None):
    gg = dbfile[gene_id]
    gn = Gene(gene_id, gg.attrs['name'], gg.attrs['chromosome'], gg.attrs['strand'],
              gg.attrs['start'], gg.attrs['end'], retrieve=True)

    majiq_config.gene_tlb[gene_id] = gn
    junction_list = {}

    num_exp = 1 if not all_exp else majiq_config.num_experiments

    for ex_grp_id in gg['exons']:
        ex_grp = gg['exons/%s' % ex_grp_id]
        ex = Exon(ex_grp.attrs['start'], ex_grp.attrs['end'], gn,
                  annot=True, isintron=False, indata=ex_grp.attrs['in_data'], retrieve=True)
        gn.exons.append(ex)
        try:
            ex.set_pcr_score(ex_grp.attrs['pcr_name'], ex_grp.attrs['score'], ex_grp.attrs['candidate'])
        except KeyError:
            pass
        if majiq_config.gcnorm:
            ex.set_gc_content_val(ex_grp.attrs['gc_content'])

        for ex_tx_id in ex_grp['tx']:
            ex_tx = ex_grp['tx/%s' % ex_tx_id]
#            TODO: Do we need trasncript? for now is None

            transcript_id = None
            ext = ExonTx(ex_tx.attrs['start'], ex_tx.attrs['end'], transcript_id, intron=False)
            ext.gene_name = gene_id
            ex.add_exon_tx(ext)
            for jj_grp_id in ex_tx["p3_junc"]:
                jj_grp = ex_tx["p3_junc/%s" % jj_grp_id]
                try:
                    junc = junction_list[jj_grp.attrs['start'], jj_grp.attrs['end']]
                except KeyError:
                    junc = Junction(jj_grp.attrs['start'], jj_grp.attrs['end'], None, None,
                                    gn, annotated=True, retrieve=True, num_exp=num_exp)
                    junc.donor_id = jj_grp.attrs['donor_id']
                    junc.acceptor_id = jj_grp.attrs['acceptor_id']
                    junction_list[jj_grp.attrs['start'], jj_grp.attrs['end']] = junc

                ext.add_3prime_junc(junc)

            for jj_grp_id in ex_tx["p5_junc"]:
                jj_grp = ex_tx["p5_junc/%s" % jj_grp_id]
                try:
                    junc = junction_list[jj_grp.attrs['start'], jj_grp.attrs['end']]
                except KeyError:
                    junc = Junction(jj_grp.attrs['start'], jj_grp.attrs['end'], None, None,
                                    gn, annotated=True, retrieve=True, num_exp=num_exp)
                    junc.donor_id = jj_grp.attrs['donor_id']
                    junc.acceptor_id = jj_grp.attrs['acceptor_id']
                    junction_list[jj_grp.attrs['start'], jj_grp.attrs['end']] = junc

                ext.add_5prime_junc(junc)

    return gn


def extract_junctions_hdf5(gene_obj, jj_grp, junction_list, annotated=True, all_exp=False):
    num_exp = 1 if not all_exp else majiq_config.num_experiments

    try:
        junc = junction_list[jj_grp.attrs['start'], jj_grp.attrs['end']]
    except KeyError:
        junc = Junction(jj_grp.attrs['start'], jj_grp.attrs['end'], None, None,
                        gene_obj, annotated=annotated, retrieve=True, num_exp=num_exp)
        junction_list[jj_grp.attrs['start'], jj_grp.attrs['end']] = junc

    return junc
