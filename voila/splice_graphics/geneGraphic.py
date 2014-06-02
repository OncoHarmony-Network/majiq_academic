import json

__author__ = 'abarrera'


class GeneGraphic(object):

    def __init__(self, name, strand, exons, junctions):
        self.name = name
        self.strand = strand
        self.exons = exons
        self.junctions = junctions

    def get_name(self):
        return self.name

    def get_strand(self):
        return self.strand

    def get_exons(self):
        return self.exons

    def get_junctions(self):
        return self.junctions

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)