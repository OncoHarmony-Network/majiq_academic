__author__ = 'abarrera'


from collections import defaultdict
import json
import numpy as np

def padding(accumulated, converted_confidence, bins, position, direction):
    """
    When there is only one way to look for the coverage.

    @param accumulated:
    @param converted_confidence:
    @param bins:
    @param position:
    @param direction:
    @return: position of the bin where the confidence is reached
    """
    bins_len = bins.size
    while accumulated < converted_confidence:
        accumulated += bins.item(position)
        position += direction
        if position < 0:
            return 0
        if position >= bins_len:
            return bins_len - 1
    return position


def find_confidence_interval(bins, meanX, confidence=0.95):
    """
    This method returns a confidence interval around the mean
    @param bins:
    @param mean:
    """

    pos_mean = max(int(meanX) - 1, 0)
    accumulated = bins.item(pos_mean)
    if accumulated >= confidence:  # TODO: what if the mean has more than the conf. interval required?
        return [pos_mean, pos_mean + 1]

    left = pos_mean - 1
    right = pos_mean + 1

    while True:
        if left < 0:  # No more variance to catch on the left
            return [0, padding(accumulated, confidence, bins, right, 1)]

        if right >= bins.size:  # No more variance to catch on the right
            return [padding(accumulated, confidence, bins, left, -1), bins.size - 1]

        # Choose the side with bigger variance
        if bins.item(left) > bins.item(right):
            accumulated += bins.item(left)
            offset = (-1, 0)
        else:
            accumulated += bins.item(right)
            offset = (0, 1)

        if accumulated >= confidence:
            return [left, right]

        left, right = left + offset[0], right + offset[1]


def find_quartiles(bins):
    """
    Iterate the bins finding how the PSI values are distributed.
    @param bins:
    @return:
    """
    accumulated = 0
    bin_index = 0
    quartiles_set = (.10, .25, .50, .75, .90)
    quartiles_values = []
    quartiles_index = 0
    while quartiles_index < len(quartiles_set) and accumulated < 1:
        accumulated += bins[bin_index]
        while quartiles_index < len(quartiles_set) and accumulated >= quartiles_set[quartiles_index]:
            quartiles_values.append(bin_index)
            quartiles_index += 1
        bin_index += 1

    return quartiles_values


def get_expected(bins):
    bins = np.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * np.arange(step / 2, 1, step)
    return np.sum(projection_prod)


def get_variance(bins, mean):
    """Compute the variance = E[X^2] - (E[X])^2"""
    bins = np.array(bins)
    step_bins = 1.0 / bins.size
    projection_prod = bins * np.arange(step_bins / 2, 1, step_bins)**2
    return np.sum(projection_prod) - mean**2


def create_array_bins(bins, confidence):
    """Compute mean, confidence intervals, quartiles and variance for a given list of bins"""
    bins = np.array(bins)
    step = 1.0 / bins.size
    mean = get_expected(bins)
    conf_interval = find_confidence_interval(bins, mean / step, confidence)
    quartiles_set = find_quartiles(bins)
    variance = get_variance(bins, mean)
    return mean, conf_interval, quartiles_set, variance


def find_excl_incl_percentages(bins, threshold):
    """
    Calculate the percentage of inclusion/exclusion given the differential bins set

    @param bins: array of bins where the sum of all elements is equal to 1
    @param threshold: the absolute value of the minimum differential delta PSI (e.g. 0.2)
    @return array of exclusion and inclusion percentages.
    """
    edges = np.linspace(-1, 1, num=len(bins))
    edges_bins = edges * bins
    bins_per_threshold = int(len(bins) * (threshold / 2))
    return [-sum(edges_bins[:int(len(bins) / 2) - bins_per_threshold]),
            sum(edges_bins[int(len(bins) / 2) + bins_per_threshold:])]


def collapse_matrix(matrix):
    """Collapse the diagonals probabilities in 1-D and return them"""
    collapse = []
    #FOR TEST matrix = array([[0, 1, 2, 3, 4, 500], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], [100, 1, 2, 3, 4, 5], ])

    matrix_corner = matrix.shape[0]+1
    for i in xrange(-matrix_corner, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)


def extract_bins_info(lsv_bins, threshold, include_lsv):
    expected_psis_bins = []
    excl_inc_perc_list = []
    collapsed_matrices = []

    for junc_matrix in lsv_bins:
        collapsed_matrices.append(collapse_matrix(np.array(junc_matrix)))

    if len(collapsed_matrices)<2:
        collapsed_matrices.append(collapsed_matrices[-1][::-1])

    for bins in collapsed_matrices:
        expected_psis_bins.append(list(bins))
        excl_inc_tuple = find_excl_incl_percentages(bins, threshold)
        excl_inc_perc_list.append(excl_inc_tuple)

        # If the delta is significant (over the threshold) or 'show-all' option, include LSV
        include_lsv = include_lsv or np.any(np.array(excl_inc_tuple)[np.array(excl_inc_tuple)>threshold])
    return expected_psis_bins, excl_inc_perc_list, include_lsv


class OrphanJunctionException(Exception):
    def __init__(self, m):
        self.message = m


class VoilaLsv(object):
    """LSV information unit managed by Voila"""
    def __init__(self, bins_list, lsv_graphic, psi1=None, psi2=None, logger=None):
        self.bins = bins_list
        self.lsv_graphic = lsv_graphic
        self.psi1 = np.array(psi1).tolist()
        self.psi2 = np.array(psi2).tolist()

        self.means = []
        self.conf_interval = []
        self.quartiles = []
        self.variances = []

        # Bins info
        self.excl_incl = None

        # Contextual info
        if lsv_graphic:
            try:
                self.set_gff3(lsv_graphic)
            except OrphanJunctionException, e:
                if logger:
                    logger.warning(e.message)
                else:
                    print "[WARNING] :: %s" % e.message
            # For LSV filtering
            self.init_categories()
        self.psi_junction = 0

    def set_bins_info(self, bins, confidence=0.95):
        self.bins = bins
        for lsv_bins in self.bins:
            m, c, q, v = create_array_bins(lsv_bins, confidence)
            self.means.append(m)
            self.conf_interval.append(c)
            self.quartiles.append(q)
            self.variances.append(v)
        return self

    def is_delta_psi(self):
        return not (self.psi2 is None)
    
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

    def get_quartiles(self):
        return self.quartiles

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

    def init_categories(self):
        # New type example: s|1e1.3o4|i|1e2.1o1
        self.categories = defaultdict()
        juns = self.get_type().split('|')
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
        self.categories['ES'] = len(exons.keys()) > 1
        self.categories['prime5'] = len(ssites) > 1
        self.categories['prime3'] = max([len(exons[e]) for e in exons]) > 1
        self.categories['njuncs'] = np.sum(['e0' not in junc for junc in juns[1:]])
        self.categories['nexons'] = len(exons.keys()) + 1
        self.categories['source'] = juns[0] == 's'
        self.categories['target'] = juns[0] == 't'
        self.categories['ir'] = ir
        if juns[0] == 't':
            self.categories['prime5'], self.categories['prime3'] = self.categories['prime3'], self.categories['prime5']

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

    def sort_bins(self, strand):
        if len(self.bins) > 2: return  #Only for 2-way LSVs
        if strand == '+' and self.get_type().startswith('t') or strand == '-' and self.get_type().startswith('t'):
            self.bins[0], self.bins[1] = self.bins[1], self.bins[0]
            self.psi_junction = 1

    def set_bed12(self, bed12_str):
        self.bed12_str = bed12_str

    def get_bed12(self):
        return self.bed12_str

    def get_gff3(self):
        return self.gff3_str

    def set_gff3(self, geneG):

        def find_exon_a5(lexonG, jidx):
            for eG in lexonG:
                if jidx in eG.get_a5_list():
                    return eG
            raise OrphanJunctionException("Orphan junction %s in lsv %s." % (repr(geneG.get_junctions()[jidx].get_coords()), self.get_id()))

        def find_exon_a3(lexonG, jidx):
            for eG in lexonG:
                if jidx in eG.get_a3_list():
                    return eG
            raise OrphanJunctionException("Orphan junction %s in lsv %s." % (repr(geneG.get_junctions()[jidx].get_coords()), self.get_id()))

        # fields = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        self.gff3_str = ''
        trans = []
        lexons = geneG.get_exons()

        chrom = geneG.chrom
        strand = geneG.strand
        gStart = lexons[0].coords[0]
        gEnd = lexons[-1].coords[1]
        first = 0
        last = 1
        if strand == '-':
            gStart = lexons[-1].coords[1]
            gEnd = lexons[0].coords[0]
            first, last = last, first

        gene_str = '\t'.join([chrom, 'majiq', 'gene', str(gStart), str(gEnd), '.', strand, '0',
                              'Name=%s;ID=%s' % (self.get_id(), self.get_id())])

        trans.append(gene_str)
        for jid, junc in enumerate(geneG.get_junctions()):
            mrna = '%s\tmajiq\tmRNA\t' % chrom
            mrna_id = '%s.%d' % (self.get_id(), jid)
            ex1 = '%s\tmajiq\texon\t' % chrom
            ex2 = '%s\tmajiq\texon\t' % chrom

            ex1G = find_exon_a5(lexons, jid)
            ex2G = find_exon_a3(lexons, jid)

            if strand == '-':
                ex1G, ex2G = ex2G, ex1G

            mrna += '%d\t%d\t' % (ex1G.get_coords()[first], ex2G.get_coords()[last])

            if self.get_type().startswith('t'):
                ex1G, ex2G = ex2G, ex1G
                ex1 += '%d\t%d\t' % (junc.get_coords()[last], ex1G.get_coords()[last])
                ex2 += '%d\t%d\t' % (ex2G.get_coords()[first], junc.get_coords()[first])
            else:
                ex1 += '%d\t%d\t' % (ex1G.get_coords()[first], junc.get_coords()[first])
                ex2 += '%d\t%d\t' % (junc.get_coords()[last], ex2G.get_coords()[last])
            mrna += '.\t%s\t0\tName=%s;Parent=%s;ID=%s' % (strand, mrna_id, self.get_id(), mrna_id)
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

        self.gff3_str = lsv_gtf

    def to_JSON(self, encoder=json.JSONEncoder):
        self.bins = np.array(self.bins).tolist()
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)
