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


def get_mean_step(bins): # TODO: rename
    bins = np.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * np.arange(step / 2, 1, step)
    return step, np.sum(projection_prod)


def get_variance(bins, mean):
    """Compute the variance = E[X^2] - (E[X])^2"""
    bins = np.array(bins)
    step_bins = 1.0 / bins.size
    projection_prod = bins * np.arange(step_bins / 2, 1, step_bins)**2
    return np.sum(projection_prod) - mean**2


def create_array_bins(bins, confidence):
    """
    Recaps bins info from data previously generated and stored in a Pickle file
    @param event_id: to access the bins associated with the event
    @param num_bins: ONLY in DEBUG (it should be retrieved from file)
    @return: a tuple with:
     *.- the mean,
     *.- the coordinates of the confidence interval [coord1, coord2].

    """
    bins = np.array(bins)
    step, mean = get_mean_step(bins)
    conf_interval = find_confidence_interval(bins, mean / step, confidence)
    quartiles_set = find_quartiles(bins)
    variance = get_variance(bins, mean)
    return mean, conf_interval, quartiles_set, variance

class Lsv(object):

    def __init__(self, lsvs_bins, lsv_meta, confidence=.2):

        means_psi_list = []
        conf_interval_list = []
        quartile_list = []
        variance_list = []

        for lsv_bins in lsvs_bins:
            m, c, q, v = create_array_bins(lsv_bins, confidence)
            means_psi_list.append(m)
            conf_interval_list.append(c)
            quartile_list.append(q)
            variance_list.append(v)

        self.id = lsv_meta[1]
        self.coords = lsv_meta[0]
        self.type = lsv_meta[2]

        # Gene info
        self.chorm = (lsv_meta[4].chrom)
        self.strand = (lsv_meta[4].strand)

        # Bins info
        self.bins = np.array(lsvs_bins).tolist()
        self.variances = variance_list
        self.means = means_psi_list
        self.conf_interval = conf_interval_list
        self.quartiles = quartile_list
        self.excl_incl = None

        # Contextual info
        self.set_extension(lsv_meta[4])
        self.set_gff3(lsv_meta[4])

        # For LSV filtering
        self.init_categories()
        self.psi_junction = 0


    def get_chrom(self):
        return self.chrom

    def get_strand(self):
        return self.strand

    def set_chrom(self, c):
        self.chrom = c

    def set_genome(self, g):
        self.genome = g

    def set_strand(self, s):
        self.strand = s

    def set_id(self, idp):
        self.id = idp

    def get_id(self):
        return self.id

    def set_type(self, t):
        self.type = t

    def get_type(self):
        return self.type

    def get_bins(self):
        return self.bins

    def get_means(self):
        return self.means

    def get_quartiles(self):
        return self.quartiles

    def set_coords(self, coords):
        self.coords = coords

    def get_coords(self):
        return self.coords

    def get_variances(self):
        return self.variances

    def set_excl_incl(self, excl_incl_set):
        self.excl_incl = excl_incl_set

    def get_excl_incl(self):
        return self.excl_incl

    def get_extension(self):
        return self.extension

    # def set_extension(self, geneG, lsv_type):
    #     def _find_lsv(coords, exon_list):
    #         for e in exon_list:
    #             if e.coords == coords:
    #                 return e
    #
    #     lsv_exon = _find_lsv(self.coords, geneG.get_exons())
    #
    #     if lsv_type.startswith('s'):
    #         self.extension = [self.coords[0], geneG.get_junctions()[lsv_exon.get_a5_list()[-1]].get_coords()[1] + 100]
    #     else:
    #         self.extension = [geneG.get_junctions()[lsv_exon.get_a3_list()[0]].get_coords()[0] - 100, self.coords[1]]

    def set_extension(self, geneG, lsv_type=None):

        if not lsv_type:
            lsv_type = self.type

        if lsv_type.startswith('s'):
            self.extension = [self.coords[0], geneG.get_exons()[-1].get_coords()[1]]
        elif lsv_type.startswith('t'):
            self.extension = [geneG.get_exons()[0].get_coords()[0], self.coords[1]]
        else:
            print "[ERROR] :: LSV type not recognized: %s" % lsv_type


    def init_categories(self):
        self.categories = defaultdict()
        j = self.type.split('|')
        ssites = set(int(s[0]) for s in j[1:])
        exons = defaultdict(list)

        for s in j[1:]:
            exs = s[1:].split('.')
            try:
                ssite = int(exs[1])
                exons[exs[0]].append(ssite)
            except IndexError:
                pass
        self.categories['ES'] = len(exons.keys()) > 1
        self.categories['prime5'] = len(ssites) > 1
        self.categories['prime3'] = max([len(exons[e]) for e in exons]) > 1
        self.categories['njuncs'] = np.sum(['e0' not in junc for junc in j[1:]])
        self.categories['nexons'] = len(exons.keys()) + 1
        self.categories['source'] = j[0] == 's'
        self.categories['target'] = j[0] == 't'

        if j[0] == 't':
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
        if strand == '+' and self.type.startswith('t') or strand == '-' and self.type.startswith('t'):
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
            print "[WARNING] :: Orphan junction %s in lsv %s." % (repr(geneG.get_junctions()[jidx]), self.id)

        def find_exon_a3(lexonG, jidx):
            for eG in lexonG:
                if jidx in eG.get_a3_list():
                    return eG
            print "[WARNING] :: Orphan junction %s in lsv %s." % (repr(geneG.get_junctions()[jidx]), self.id)

        # fields = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
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
                              'Name=%s;ID=%s' % (self.id, self.id)])

        trans.append(gene_str)
        for jid, junc in enumerate(geneG.get_junctions()):
            mrna = '%s\tmajiq\tmRNA\t' % chrom
            mrna_id = '%s.%d' % (self.id, jid)
            ex1 = '%s\tmajiq\texon\t' % chrom
            ex2 = '%s\tmajiq\texon\t' % chrom

            ex1G = find_exon_a5(lexons, jid)
            ex2G = find_exon_a3(lexons, jid)

            if strand == '-':
                ex1G, ex2G = ex2G, ex1G

            mrna += '%d\t%d\t' % (ex1G.get_coords()[first], ex2G.get_coords()[last])

            if self.type.startswith('t'):
                ex1G, ex2G = ex2G, ex1G
                ex1 += '%d\t%d\t' % (junc.get_coords()[last], ex1G.get_coords()[last])
                ex2 += '%d\t%d\t' % (ex2G.get_coords()[first], junc.get_coords()[first])
            else:
                ex1 += '%d\t%d\t' % (ex1G.get_coords()[first], junc.get_coords()[first])
                ex2 += '%d\t%d\t' % (junc.get_coords()[last], ex2G.get_coords()[last])
            mrna += '.\t%s\t0\tName=%s;Parent=%s;ID=%s' % (strand, mrna_id, self.id, mrna_id)
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
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)