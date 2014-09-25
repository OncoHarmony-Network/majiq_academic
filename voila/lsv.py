from collections import defaultdict
import json
import numpy
from splice_graphics.geneGraphic import GeneGraphic


class Lsv(object):

    def __init__(self, input_data=None):

        self.id = input_data['id']
        self.name = input_data['name']
        self.type = input_data['type'] # TODO: COnsider changing type to LSV type
        # self.bins = input_data['bins_list']
        self.bins = numpy.array(input_data['bins_list']).tolist()

        self.coords = input_data['coords']
        self.matrix_link = input_data['matrix_link']
        self.variances = input_data['variance_list']
        self.means = input_data['mean_psi']
        self.conf_interval = input_data['conf_interval']
        self.quartiles = input_data['quartiles']
        self.excl_incl = None
        self.extension = None

        # For LSV filtering
        self.set_categories()

    def get_id(self):
        return self.id

    def get_name(self):
        return self.name

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

    def get_matrix_link(self):
        return self.matrix_link

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

    def set_extension(self, geneG, lsv_type):

        if lsv_type.startswith('s'):
            self.extension = [self.coords[0], geneG.get_exons()[-1].get_coords()[1]]
        elif lsv_type.startswith('t'):
            self.extension = [geneG.get_exons()[0].get_coords()[0], self.coords[1]]
        else:
            print "ERRORR! LSV type not recognized: %s" % lsv_type


    def set_categories(self):
        self.categories = defaultdict()
        j = self.type.split('|')
        ssites = set(int(s[0]) for s in j[1:])
        exons = defaultdict(list)
        for s in j[1:]:
            exs = s[1:].split('.')
            try:
                exons[exs[0]].append(int(exs[1]))
            except IndexError:
                pass
        self.categories['ES'] = len(exons.keys()) > 1
        self.categories['prime5'] = len(ssites) > 1
        self.categories['prime3'] = max([len(exons[e]) for e in exons]) > 1
        self.categories['njuncs'] = len(j[1:])
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

    def is_three_prime(self):
        pass

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)