import json

__author__ = 'abarrera'


class GeneGraphic(object):

    def __init__(self, name, strand, exons, junctions, chrom=None):
        self.name = name
        self.strand = strand
        self.exons = exons
        self.junctions = junctions
        self.chrom = chrom

    def get_name(self):
        return self.name

    def get_strand(self):
        return self.strand

    def get_exons(self):
        return self.exons

    def get_junctions(self):
        return self.junctions

    def get_chrom(self):
        return self.chrom

    def get_coords(self):
        if self.strand == '+':
            return [self.exons[0].coords[0], self.exons[-1].coords[1]]
        if self.strand == '-':
            return [self.exons[-1].coords[1], self.exons[0].coords[0]]

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)