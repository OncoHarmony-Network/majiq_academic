import json


__author__ = 'abarrera'


class ExonGraphic(object):

    def __init__(self, a3, a5, coords, type_exon, coords_extra=[], intron_retention=None, lsv_type=0):
        self.a3 = a3
        self.a5 = a5
        self.coords = coords
        self.type_exon = type_exon
        self.coords_extra = coords_extra
        self.intron_retention = intron_retention
        self.lsv_type = lsv_type

    def get_a3_list(self):
        return self.a3

    def get_a5_list(self):
        return self.a5

    def get_coords(self):
        return self.coords

    def get_type(self):
        return self.type_exon

    def get_coords_extra(self):
        return self.coords_extra

    def get_intron_retention(self):
        return self.intron_retention

    def get_lsv_type(self):
        return self.lsv_type

    def set_intron_retention(self, intron_retention):
        self.intron_retention = intron_retention

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)