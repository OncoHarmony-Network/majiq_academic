import json

__author__ = 'abarrera'


class JunctionGraphic(object):
    def __init__(self, coords, type_junction, num_reads):
        self.coords = coords
        self.type_junction = type_junction
        self.num_reads = num_reads

    def get_coords(self):
        return self.coords

    def get_type(self):
        return self.type_junction

    def get_num_reads(self):
        return self.num_reads

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)