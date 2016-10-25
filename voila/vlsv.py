from collections import defaultdict
import json
import numpy as np


def get_expected_dpsi(bins):
    return sum(np.array(bins) * np.arange(-1+1./len(bins), 1., 2./len(bins)))


def get_expected_psi(bins):
    bins = np.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * np.arange(step / 2, 1, step)
    return np.sum(projection_prod)


def collapse_matrix(matrix):
    """Collapse the diagonals probabilities in 1-D and return them"""
    collapse = []
    matrix_corner = matrix.shape[0]
    for i in xrange(-matrix_corner+1, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)


def _find_delta_border(V, numbins):
    """Finds the border index to which a V corresponds in its delta_space given the number of bins the matrix will have"""
    delta_space = list(np.linspace(-1, 1, num=numbins + 1))
    delta_space.pop(0)  # first border to the left is -1, and we are not interested in it
    # get the index position that corresponds to the V threshold
    for i, value in enumerate(delta_space):
        if value > V:
            return i
    # if nothing hit, V = 1
    return numbins


def matrix_area(matrix, V=0.2, absolute=True, collapsed_mat=False):
    """Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute"""
    collapse = matrix
    if not collapsed_mat:
        collapse = collapse_matrix(matrix)
    # get the delta psi histogram borders based on the size of 'collapse'
    border = _find_delta_border(V, collapse.shape[0])
    # grab the values inside the area of interest
    area = []
    if V < 0:
        area.append(collapse[0:border + 1])
        if absolute:  # if absolute V, pick the other side of the array
            area.append(collapse[-border - 1:])
    else:
        area.append(collapse[border:])
        if absolute and border != 0:  # if absolute V, pick the other side of the array
            area.append(collapse[0:len(collapse) - border])
    return sum(area)


class OrphanJunctionException(Exception):
    def __init__(self, m):
        self.message = m


class VoilaLsv(object):
    """LSV information unit managed by Voila"""

    @classmethod
    def to_gff3(cls, vlsv):

        def find_exon_a5(lexonG, jidx):
            for eG in lexonG:
                if jidx in eG.get_a5_list():
                    return eG
            raise OrphanJunctionException("Orphan junction %s in lsv %s." % (repr(vlsv.get_lsv_graphic().get_junctions()[jidx].get_coords()), lsvId))

        def find_exon_a3(lexonG, jidx):
            for eG in lexonG:
                if jidx in eG.get_a3_list():
                    return eG
            raise OrphanJunctionException("Orphan junction %s in lsv %s." % (repr(vlsv.get_lsv_graphic().get_junctions()[jidx].get_coords()), lsvId))

        # fields = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        trans = []
        lexons = vlsv.get_lsv_graphic().get_exons()
        lsvId = vlsv.get_id()
        chrom = vlsv.get_chrom()
        strand = vlsv.get_strand()
        gStart = lexons[0].coords[0]
        gEnd = lexons[-1].coords[1]
        first = 0
        last = 1
        if strand == '-':
            gStart = lexons[-1].coords[1]
            gEnd = lexons[0].coords[0]
            first, last = last, first

        gene_str = '\t'.join([chrom, 'majiq', 'gene', str(gStart), str(gEnd), '.', strand, '0',
                              'Name=%s;ID=%s' % (lsvId, lsvId)])

        trans.append(gene_str)
        for jid, junc in enumerate(vlsv.get_lsv_graphic().get_junctions()):
            mrna = '%s\tmajiq\tmRNA\t' % chrom
            mrna_id = '%s.%d' % (lsvId, jid)
            ex1 = '%s\tmajiq\texon\t' % chrom
            ex2 = '%s\tmajiq\texon\t' % chrom

            ex1G = find_exon_a5(lexons, jid)
            ex2G = find_exon_a3(lexons, jid)

            if strand == '-':
                ex1G, ex2G = ex2G, ex1G

            mrna += '%d\t%d\t' % (ex1G.get_coords()[first], ex2G.get_coords()[last])

            if vlsv.get_type().startswith('t'):
                ex1G, ex2G = ex2G, ex1G
                ex1 += '%d\t%d\t' % (junc.get_coords()[last], ex1G.get_coords()[last])
                ex2 += '%d\t%d\t' % (ex2G.get_coords()[first], junc.get_coords()[first])
            else:
                ex1 += '%d\t%d\t' % (ex1G.get_coords()[first], junc.get_coords()[first])
                ex2 += '%d\t%d\t' % (junc.get_coords()[last], ex2G.get_coords()[last])
            mrna += '.\t%s\t0\tName=%s;Parent=%s;ID=%s' % (strand, mrna_id, lsvId, mrna_id)
            ex1 += '.\t%s\t0\tName=%s.lsv;Parent=%s;ID=%s.lsv' % (strand, mrna_id, mrna_id, mrna_id)
            ex2 += '.\t%s\t0\tName=%s.ex;Parent=%s;ID=%s.ex' % (strand, mrna_id, mrna_id, mrna_id)
            trans.append(mrna)

            trans.append(ex1)
            trans.append(ex2)
        else:
            if strand == '-':
                for ii, t in enumerate(trans):
                    t_fields = t.split()
                    t_fields[3],  t_fields[4] = t_fields[4], t_fields[3]
                    trans[ii] = '\t'.join(t_fields)

            lsv_gtf = '\n'.join(trans)

        return lsv_gtf

    @classmethod
    def init_categories(cls, lsv_type):
        categories = defaultdict()
        juns = lsv_type.split('|')
        ir = 'i' in juns
        try:
            juns.remove('i')
        except ValueError:
            pass  # No sign of intron retention
        ssites = set(int(s[0]) for s in juns[1:])
        exons = defaultdict(list)

        for s in juns[1:]:
            exs = s[1:].split('.')
            try:
                ssite = int(exs[1].split('o')[0])
                exons[exs[0]].append(ssite)
            except IndexError:
                pass
        categories['ES'] = len(exons.keys()) > 1
        categories['prime5'] = len(ssites) > 1
        categories['prime3'] = max([len(exons[e]) for e in exons]) > 1
        categories['njuncs'] = np.sum(['e0' not in junc for junc in juns[1:]])
        categories['nexons'] = len(exons.keys()) + 1
        categories['source'] = juns[0] == 's'
        categories['target'] = juns[0] == 't'
        categories['ir'] = ir
        if juns[0] == 't':
            categories['prime5'], categories['prime3'] = categories['prime3'], categories['prime5']
        return categories

    def __init__(self, bins_list, lsv_graphic, psi1=None, psi2=None, logger=None):
        self.lsv_graphic = lsv_graphic
        self.psi1 = psi1
        self.psi2 = psi2

        self.means = []
        self.variances = []

        self.bins = bins_list
        self.excl_incl = []

        if self.is_delta_psi():  # Store collapsed matrix to save some space
            self.bins = [collapse_matrix(np.array(lsv_bins)) for lsv_bins in self.bins]

        if len(self.bins) == 1:
            self.bins.append(self.bins[-1][::-1])  # Recreate complementary junction in binary LSV

        for lsv_bins in self.bins:
            if self.is_delta_psi():
                self.means.append(get_expected_dpsi(lsv_bins))
                if self.means[-1] < 0:
                    self.excl_incl.append([-self.means[-1], 0])
                else:
                    self.excl_incl.append([0, self.means[-1]])
            else:
                self.means.append(get_expected_psi(np.array(lsv_bins)))
                step_bins = 1.0 / lsv_bins.size
                projection_prod = lsv_bins * np.arange(step_bins / 2, 1, step_bins)**2
                self.variances.append(np.sum(projection_prod) - self.means[-1]**2)

        # For LSV filtering
        if lsv_graphic:
            self.categories = VoilaLsv.init_categories(self.get_type())
        self.psi_junction = 0

    def is_delta_psi(self):
        return not (self.psi2 is None)

    def get_lsv_graphic(self):
        return self.lsv_graphic

    def get_chrom(self):
        return self.lsv_graphic.chrom

    def get_strand(self):
        return self.lsv_graphic.strand

    def set_chrom(self, c):
        self.lsv_graphic.chrom = c

    def set_strand(self, s):
        self.lsv_graphic.strand = s

    def set_id(self, idp):
        self.lsv_graphic.id = idp

    def get_id(self):
        return self.lsv_graphic.id

    def get_gene_name(self):
        return self.lsv_graphic.name

    def set_type(self, t):
        self.lsv_graphic.type = t

    def get_type(self):
        return self.lsv_graphic.type

    def get_bins(self):
        return self.bins

    def set_means(self, m):
        self.means = m

    def get_means(self):
        return self.means

    def set_coords(self, coords):
        self.lsv_graphic.coords = coords

    def get_coords(self):
        return self.lsv_graphic.coords

    def get_variances(self):
        return self.variances

    def set_excl_incl(self, excl_incl_set):
        self.excl_incl = excl_incl_set

    def get_excl_incl(self):
        return self.excl_incl

    def get_extension(self):
        return [self.lsv_graphic.get_exons()[0].get_coords()[0], self.lsv_graphic.get_exons()[-1].get_coords()[1]]

    def get_categories(self):
        return self.categories

    def njuncs(self):
        return self.categories['njuncs']

    def nexons(self):
        return self.categories['nexons']

    def categories2css(self):
        css_cats = []
        for c in self.categories:
            if type(self.categories[c]) == bool and self.categories[c]:
                css_cats.append(c)
        return ' '.join(css_cats)

    def set_bed12(self, bed12_str):
        self.bed12_str = bed12_str

    def get_bed12(self):
        return self.bed12_str

    def get_gff3(self, logger=None):
        try:
            return VoilaLsv.to_gff3(self)
        except OrphanJunctionException, e:
            if logger:
                logger.warning(e.message)
            else:
                print "[WARNING] :: %s" % e.message

    def to_JSON(self, encoder=json.JSONEncoder):
        self.bins = np.array(self.bins).tolist()
        if self.is_delta_psi():
            self.psi1 = np.array(self.psi1).tolist()
            self.psi2 = np.array(self.psi2).tolist()
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def is_lsv_changing(self, thres):
        means = np.array(self.get_means())
        # TODO: should we check that pos and neg are kind of matched?
        return max(means[means > 0].sum(), abs(means[means < 0].sum())) >= thres

