import pickle
from majiq.src.config import Config
from voila import constants as voila_const
from voila.splice_graphics import ExonGraphic, LsvGraphic, JunctionGraphic
from majiq.src.constants import *
import majiq.src.normalize as majiq_norm
import collections
from majiq.src.sample import sample_from_junctions

__author__ = 'jordi@biociphers.org'

SSOURCE = 'source'
STARGET = 'target'

quant_lsv = collections.namedtuple('quant_lsv', 'id type coverage')


def dump_bin_file(data, filename):
    with open(filename, 'wb') as ofp:
        fast_pickler = pickle.Pickler(ofp, protocol=2)
        # fast_pickler.fast = 1
        fast_pickler.dump(data)


class InvalidLSV(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class LSV(object):
    def __init__(self, exon, lsv_id, junctions, lsv_type):

        majiq_config = Config()
        if lsv_type != SSOURCE and lsv_type != STARGET:
            raise InvalidLSV('Incorrect LSV type %s' % lsv_type)
        self.coords = exon.get_coordinates()
        self.id = lsv_id
        junction_list = [x for x in junctions if x is not None and
                                                 x.get_donor() is not None and
                                                 x.get_acceptor() is not None]
        if majiq_config.simplify:
            jj_set = set()
            for exp_idx in range(majiq_config.num_experiments):
                cover = [float(junc.get_coverage_sum(exp_idx)) for junc in junction_list]
                if sum(cover) == 0:
                    continue
                bool_map = [majiq_config.simplify_type == SIMPLIFY_ALL or
                            (junc.is_annotated() and majiq_config.simplify_type == SIMPLIFY_DB) or
                            (not junc.is_annotated() and majiq_config.simplify_type == SIMPLIFY_DENOVO)
                            for junc in junction_list]

                jj_set = jj_set.union(set([junc for eidx, junc in enumerate(junction_list)
                                           if cover[eidx]/sum(cover) >= majiq_config.simplify_threshold and bool_map[eidx]]))
            junction_list = list(jj_set)

        n_viable_juncs = len(junction_list)
        if n_viable_juncs < 2:
            raise InvalidLSV('Not enought junctions')
        self.type = lsv_type
        self.exon = exon

        self.intron_retention = False
        for jj in junction_list:
            x1 = jj.get_acceptor()
            x2 = jj.get_donor()

            if x1 is None or x2 is None:
                continue
            if x1.is_intron() or x2.is_intron():
                self.intron_retention = True
                break
        try:
            tlb_junc = {}
            self.ext_type = self.set_type(junction_list, tlb_junc)
            if self.ext_type == 'intron':
                raise InvalidLSV('Auto junction found')
        except:
            raise InvalidLSV('Problematic Type')

        self.junctions = []
        order = self.ext_type.split('|')[1:]
        for idx, jj in enumerate(order):
            if jj[-2:] == 'e0':
                continue
            self.junctions.append(junction_list[tlb_junc[jj]])
        del junction_list
        #self.visual = list()
        # for exp_idx in range(majiq_config.num_experiments):
        #     self.visual.append(self.get_visual_lsv(self.junctions, exp_idx))
        # self.visual = np.array(self.visual)

    def get_coordinates(self):
        return self.coords

    def get_junctions_list(self):
        return self.junctions

    def is_Ssource(self):
        return bool(self.type == SSOURCE)

    def is_Starget(self):
        return bool(self.type == STARGET)

    def has_pcr_score(self):
        return not self.exon.get_pcr_score() is None

    def get_pcr_score(self):
        return self.exon.get_pcr_score()

    def get_strand(self):
        return self.exon.get_strand()

    def get_chromosome(self):
        return self.exon.get_gene().get_chromosome()

    def get_visual(self, exp_idx):
        return self.visual[exp_idx]

    def set_type(self, jlist, tlb_junc):
        ex_id = self.exon.get_id()
        strand = self.get_strand()
        rev = (strand == '-')

        if self.type == SSOURCE:
            spsite = sorted(set(self.exon.ss_5p_list), reverse=rev)
        else:
            spsite = sorted(set(self.exon.ss_3p_list), reverse=rev)
        ex_set = set()
        skip = False
        for junc in jlist:
            jdonor = junc.get_donor()
            jacceptor = junc.get_acceptor()
            if self.type == SSOURCE:
                lsv_exon = jdonor
                if lsv_exon.get_id() != ex_id:
                    skip = True
                    break

                if not jacceptor is None:
                    if not jacceptor.is_intron():
                        ex_set.add(jacceptor.get_id())
            else:
                lsv_exon = jacceptor
                if jacceptor is None:
                    continue
                if lsv_exon.get_id() != ex_id:
                    skip = True
                    break
                if not jdonor is None:
                    if not jdonor.is_intron():
                        ex_set.add(jdonor.get_id())

        if skip:
            return 'intron'

        ex_list = sorted(list(ex_set), reverse=rev)

        if (self.type == SSOURCE and strand == '+') or (self.type == STARGET and strand == '-'):
            ext_type = "s"
        else:
            ext_type = "t"

        type_set = set()
        for jidx, junc in enumerate(jlist):
            jdonor = junc.get_donor()
            jacceptor = junc.get_acceptor()
            if self.type == SSOURCE:
                if jacceptor is None:
                    exs3 = ''
                    ex = '0'
                    jtype = "|%se%s" % (spsite.index(junc.start) + 1, ex)
                elif jacceptor.is_intron():
                    jtype = "|i"
                else:
                    s3 = sorted(list(set(jacceptor.ss_3p_list)), reverse=rev)
                    ex1 = ex_list.index(jacceptor.get_id()) + 1
                    try:
                        ex = '%s.%so%s' % (ex1, s3.index(junc.end) + 1, len(s3))
                    except Exception as e:
                        print("ERRORRR", ex_id, e)
                        raise e
                    jtype = "|%se%s" % (spsite.index(junc.start) + 1, ex)
            else:
                if jdonor is None:
                    exs5 = ''
                    ex = '0'
                    jtype = "|%se%s" % (spsite.index(junc.end) + 1, ex)
                elif jdonor.is_intron():
                    jtype = "|i"
                else:
                    s5 = sorted(list(set(jdonor.ss_5p_list)), reverse=rev)
                    ex1 = ex_list.index(jdonor.get_id()) + 1
                    ex = '%s.%so%s' % (ex1, s5.index(junc.start) + 1, len(s5))
                    jtype = "|%se%s" % (spsite.index(junc.end) + 1, ex)
            type_set.add(jtype)
            tlb_junc[jtype[1:]] = jidx
        for tt in sorted(list(type_set)):
            ext_type += tt

        return ext_type

    def get_visual_lsv(self, junction_list, exp_idx):
        junc_list = []
        junc_l = []
        lsv_exon_list = [self.exon]
        alt_empty_ends = []
        alt_empty_starts = []

        for jj in junction_list:
            jdonor = jj.get_donor()
            jacceptor = jj.get_acceptor()
            cc = jj.get_coordinates()
            if jdonor is None:
                alt_empty_ends.append(cc[1])
                continue
            if jacceptor is None:
                alt_empty_starts.append(cc[0])
                continue

            if jacceptor != self.exon:
                lsv_exon_list.append(jacceptor)
            if jdonor != self.exon:
                lsv_exon_list.append(jdonor)

            if jj.annotated and jj.get_read_num(exp_idx) == 0:
                jtype = 2
            elif jj.annotated and jj.get_read_num(exp_idx) > 0:
                jtype = 0
            elif not jj.annotated and jj.get_read_num(exp_idx) > 0:
                jtype = 1
            else:
                jtype = 1
                # continue

            ir_type = 0
            if jj.get_donor().is_intron():
                ir_type = voila_const.IR_TYPE_START
            elif jj.get_acceptor().is_intron():
                ir_type = voila_const.IR_TYPE_END

            junc_l.append(jj.get_coordinates())
            junc_list.append(JunctionGraphic(jj.get_coordinates()[0], jj.get_coordinates()[1],
                                             junction_type_list=[], reads_list=[],
                                             transcripts=jj.get_transcript_list(), intron_retention=ir_type))
        junc_l = np.asarray(junc_l)
        lsv_exon_list.sort()
        exon_list = []
        for ex in lsv_exon_list:
            cc = ex.get_coordinates()
            a3 = []
            alt_start = []
            for ss3 in set(ex.ss_3p_list):
                if ss3 in alt_empty_starts:
                    alt_start.append(ss3)
                    # continue
                for jidx, jjl in enumerate(junc_l):
                    if ss3 == jjl[1]:
                        a3.append(jidx)

            a5 = []
            alt_ends = []
            for ss5 in set(ex.ss_5p_list):
                if ss5 in alt_empty_starts:
                    alt_ends.append(ss5)
                    # continue
                for jidx, jjl in enumerate(junc_l):
                    if ss5 == jjl[0]:
                        a5.append(jidx)

            ex_reads = ex.get_total_read_num(exp_idx)

            if ex.is_annotated() and ex_reads == 0.0:
                visual_type = 2
            elif ex.is_annotated() and ex_reads > 0.0:
                visual_type = 0
            elif not ex.is_annotated() and ex_reads > 0.0:
                visual_type = 1
            else:
                visual_type = 1
                # continue
            extra_coords = []
            if ex.is_annotated():
                if ex.start < ex.db_coord[0]:
                    extra_coords.append([ex.start, ex.db_coord[0] - 1])
                if ex.end > ex.db_coord[1]:
                    extra_coords.append([ex.db_coord[1] + 1, ex.end])

            eg = ExonGraphic(a3, a5, start=cc[0], end=cc[1], exon_type_list=[], coords_extra=extra_coords,
                             intron_retention=ex.get_ir(), alt_starts=alt_start, alt_ends=alt_ends)
            exon_list.append(eg)

        splice_lsv = LsvGraphic(lsv_type=self.ext_type, start=self.coords[0], end=self.coords[1], lsv_id=self.id,
                                name=self.exon.get_gene().get_name(), chromosome=self.get_chromosome(),
                                strand=self.get_strand(), exons=exon_list, junctions=junc_list)
        return splice_lsv

    def is_equivalent(self, variant):

        jlist1 = sorted(self.junctions)
        jlist2 = sorted(variant.junctions)
        if self.type == variant.type:
            return False
        return np.array_equal(jlist1, jlist2)

    def contained(self, variant):
        #TODO: Change for set subset
        res = True
        jlist1 = sorted(self.junctions)
        jlist2 = sorted(variant.junctions)
        if self.type == variant.type:
            return False
        if np.array_equal(jlist1, jlist2):
            if (self.get_strand() == '+' and self.type == SSOURCE) or \
                    (self.get_strand() == '-' and self.type == STARGET):
                res = False
            else:
                res = True
        else:
            for jj1 in jlist1:
                if jj1 not in jlist2:
                    res = False
                    break

        return res

#     def to_hdf5(self, hdf5grp, lsv_idx, exp_idx, fitfunc_r=1, gc_vfunc=None):
#         majiq_config = Config()
#         try:
#             njunc = len(self.junctions)
#             cover = np.zeros(shape=(njunc, (majiq_config.readLen - 16) + 1),
#                              dtype=np.float)
#
#             pvalue_limit = 0.0000001
#             for idx, junc in enumerate(self.junctions):
#                 if junc.get_index() != -1:
#                     cover[idx] = junc.get_coverage()[exp_idx]
#                     if majiq_config.gcnorm and junc.get_gc_content(exp_idx).sum() > 0:
#                         vals = gc_vfunc(junc.get_gc_content(exp_idx))
#                         cover[idx] = np.multiply(cover[idx], vals)
#                     if pvalue_limit >= 0:
#                         cover[idx] = majiq_norm.mark_stacks_per_junc(cover[idx], fitfunc_r, pvalue_limit)
#
# #            if lsv_idx + njunc > majiq_config.nrandom_junctions:
#             if lsv_idx + njunc > 2:
#                 shp = hdf5grp[JUNCTIONS_DATASET_NAME].shape
#                 shp_new = shp[0] + NRANDOM_JUNCTIONS
#                 hdf5grp[JUNCTIONS_DATASET_NAME].resize((shp_new, shp[1]))
#
#             hdf5grp[JUNCTIONS_DATASET_NAME][lsv_idx:lsv_idx+njunc] = cover
#
#             h_lsv = hdf5grp.create_group("LSVs/%s" % self.id)
#             # h_lsv.attrs['coords'] = self.coords
#             h_lsv.attrs['id'] = self.id
#             h_lsv.attrs['type'] = self.ext_type
#             # h_lsv.attrs['coverage'] = hdf5grp[JUNCTIONS_DATASET_NAME].regionref[lsv_idx:lsv_idx + njunc]
#             h_lsv.attrs['coverage_index'] = (lsv_idx, lsv_idx + njunc)
#
#             vh_lsv = h_lsv.create_group('visual')
#             self.get_visual(exp_idx).to_hdf5(vh_lsv)
#
#         except:
#             print("HDF5 ERROR", self.id, cover.shape, hdf5grp[JUNCTIONS_DATASET_NAME].shape)
#             raise
#
#         return lsv_idx + njunc

    def to_hdf5(self, hdf5grp, gc_vfunc, fitfunc_r, exp, exp_idx):

        majiq_config = Config()
        lsv_idx = hdf5grp.attrs['lsv_idx']
        njunc = len(self.junctions)
        cover = np.zeros(shape=(njunc, (majiq_config.readLen - 16) + 1), dtype=np.float)

        pvalue_limit = majiq_config.markstacks
        for idx, junc in enumerate(self.junctions):
            if junc.get_index() != -1:
                cover[idx] = exp[JUNCTIONS_DATASET_NAME][junc.idx]
                if majiq_config.gcnorm:
                    mm = exp[JUNCTIONS_GC_CONTENT][junc.idx].sum()
                    if mm > 0:
                        vals = gc_vfunc(mm)
                        cover[idx] = np.multiply(cover[idx], vals)
                if pvalue_limit >= 0:
                    cover[idx] = majiq_norm.mark_stacks_per_junc(cover[idx], fitfunc_r, pvalue_limit)

        s_lsv = sample_from_junctions(junction_list=cover,
                                      m=majiq_config.m,
                                      k=majiq_config.k,
                                      fitted_one_over_r=fitfunc_r,
                                      debug=majiq_config.debug)
        lsv_trs = np.array([cover.sum(axis=1), np.count_nonzero(cover, axis=1)]).T

        njunc = s_lsv.shape[0]
        if lsv_idx + njunc > 2:
            shp = hdf5grp[JUNCTIONS_DATASET_NAME].shape
            shp_new = shp[0] + 5000
            hdf5grp[JUNCTIONS_DATASET_NAME].resize((shp_new, shp[1]))
            hdf5grp['junc_cov'].resize((shp_new, 2))

        hdf5grp[JUNCTIONS_DATASET_NAME][lsv_idx:lsv_idx + njunc] = s_lsv
        # print('%s::%s' % (vals['junc_attr'], type(vals['junc_attr'])))
        hdf5grp['junc_cov'][lsv_idx:lsv_idx + njunc] = lsv_trs

        h_lsv = hdf5grp.create_group("LSVs/%s" % self.id)
        h_lsv.attrs['id'] = self.id
        h_lsv.attrs['type'] = self.ext_type
        h_lsv.attrs['coverage'] = [lsv_idx, lsv_idx + njunc]
        # TODO: CHECK
        vh_lsv = h_lsv.create_group('visual')
        #self.get_visual(exp_idx).to_hdf5(h_lsv)
        self.get_visual_lsv(self.junctions, exp_idx).to_hdf5(vh_lsv)

        # lsv_idx = boots_write(hdf5grp, vals, lsv_idx)
        hdf5grp.attrs['lsv_idx'] = lsv_idx + njunc
        hdf5grp.attrs['num_lsvs'] = hdf5grp.attrs['num_lsvs'] + 1

    def to_queue(self, gc_vfunc, fitfunc_r, exp, exp_idx):
        majiq_config = Config()

        njunc = len(self.junctions)
        cover = np.zeros(shape=(njunc, (majiq_config.readLen - 16) + 1), dtype=np.float)

        pvalue_limit = majiq_config.markstacks
        for idx, junc in enumerate(self.junctions):
            if junc.get_index() != -1:
                cover[idx] = exp[JUNCTIONS_DATASET_NAME][junc.idx]
                if majiq_config.gcnorm:
                    mm = exp[JUNCTIONS_GC_CONTENT][junc.idx].sum()
                    if mm > 0:
                        vals = gc_vfunc(mm)
                        cover[idx] = np.multiply(cover[idx], vals)
                if pvalue_limit >= 0:
                    cover[idx] = majiq_norm.mark_stacks_per_junc(cover[idx], fitfunc_r, pvalue_limit)

        # m_lsv, var_lsv, s_lsv=sample_from_junctions(junction_list=cover,
        s_lsv = sample_from_junctions(junction_list=cover,
                                      m=majiq_config.m,
                                      k=majiq_config.k,
                                      fitted_one_over_r=fitfunc_r,
                                      debug=majiq_config.debug)
        lsv_trs = np.array([cover.sum(axis=1), np.count_nonzero(cover, axis=1)]).T
        # print(lsv_trs, np.count_nonzero(cover, axis=1), np.count_nonzero(cover, axis=0))
        #
        # print(cover)
        vals = {'samples': s_lsv, 'id': self.id, 'type': self.ext_type, 'junc_attr': lsv_trs,
                'lsv_graphic': self.get_visual_lsv(self.junctions, exp_idx)}
        return vals


def extract_se_events(list_lsv_per_gene):
    sslist = list_lsv_per_gene[0]
    stlist = list_lsv_per_gene[1]

    for ss in sslist:
        slist = ss.junction_list
        if len(slist) != 2:
            continue
        if slist[0].acceptor == slist[1].acceptor:
            continue

        for st in stlist:
            tlist = st.junction_list
            if len(tlist) != 2:
                continue
            if tlist[0].donor == tlist[1].donor:
                continue

            for ii in range(2):
                for jj in range(2):
                    if slist[ii] == tlist[jj]:
                        sindx = 1 - ii
                        tindx = 1 - jj
                        break
                else:
                    continue
                break
            else:
                continue
            if slist[sindx].acceptor == tlist[tindx].donor:
                C2 = st.coords
                break
        else:
            continue

        C1 = ss.coords
        A = slist[sindx].acceptor.get_coordinates()
        # ret_list.append( (C1.)


def extract_gff(list_lsv, out_dir):
    gtf = set()
    for name, lsv_l in list_lsv.items():
        for lsv in lsv_l:
            trans = []
            jlist = lsv.junctions
            lsv_coord = lsv.get_coordinates()

            gne = jlist[0].get_gene()
            chrom = gne.get_chromosome()
            strand = gne.get_strand()
            gene = '%s\tscript\tgene\t' % chrom
            if lsv.type == SSOURCE:
                if jlist[-1].get_acceptor() is None:
                    continue
                gene += '%d\t%d\t' % (lsv_coord[0], jlist[-1].get_acceptor().get_coordinates()[1])
            else:
                if jlist[0].get_donor() is None:
                    continue
                gene += '%d\t%d\t' % (jlist[0].get_donor().get_coordinates()[0], lsv_coord[1])

            gene += '.\t%s\t.\tName=%s;Parent=%s;ID=%s' % (strand, lsv.id, lsv.id, lsv.id)
            trans.append(gene)
            for jidx, junc in enumerate(jlist):
                mrna = '%s\tscript\tmRNA\t' % chrom
                mrna_id = '%s.%d' % (lsv.id, jidx)
                ex1 = '%s\tscript\texon\t' % chrom
                ex2 = '%s\tscript\texon\t' % chrom
                if lsv.type == SSOURCE:
                    if junc.get_acceptor() is None:
                        break
                    excoord = junc.get_acceptor().get_coordinates()
                    variant = junc.get_coordinates()
                    mrna += '%d\t%d\t' % (lsv_coord[0], excoord[1])
                    ex1 += '%d\t%d\t' % (lsv_coord[0], variant[0])
                    ex2 += '%d\t%d\t' % (variant[1], excoord[1])
                else:
                    if junc.get_donor() is None:
                        break
                    excoord = junc.get_donor().get_coordinates()
                    variant = junc.get_coordinates()
                    mrna += '%d\t%d\t' % (excoord[0], lsv_coord[1])
                    ex1 += '%d\t%d\t' % (variant[1], lsv_coord[1])
                    ex2 += '%d\t%d\t' % (excoord[0], variant[0])
                mrna += '.\t%s\t.\tName=%s;Parent=%s;ID=%s' % (strand, mrna_id, lsv.id, mrna_id)
                ex1 += '.\t%s\t.\tName=%s.lsv;Parent=%s;ID=%s.lsv' % (strand, mrna_id, mrna_id, mrna_id)
                ex2 += '.\t%s\t.\tName=%s.ex;Parent=%s;ID=%s.ex' % (strand, mrna_id, mrna_id, mrna_id)
                trans.append(mrna)
                trans.append(ex1)
                trans.append(ex2)
            else:
                lsv_gtf = '\n'.join(trans)
                gtf.add(lsv_gtf)

    gtf = sorted(gtf)
    fname = '%s/temp_gff.pkl' % out_dir
    with open(fname, 'w+b') as ofp:
        dump_bin_file(gtf, fname)

    return gtf


def new_lsv_definition(exon, jlist, lsv_type):
    coords = exon.get_coordinates()
    lsv_id = "%s:%d-%d:%s" % (exon.get_gene().get_id(), coords[0], coords[1], lsv_type)
    return LSV(exon, lsv_id, jlist, lsv_type)


def print_lsv_extype(list_lsv, filename):
    fp = open(filename, 'w+')
    print(list_lsv.shape)
    for idx in range(list_lsv.shape[0]):
        lsv = list_lsv[idx]
        fp.write("%s\n" % lsv.type)
    fp.close()