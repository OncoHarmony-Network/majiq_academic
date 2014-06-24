import json


class Lsv(object):

    def __init__(self, input_data=None):

        self.id = input_data['id']
        self.name = input_data['name']
        self.type = input_data['type']
        self.bins = input_data['bins_list'].tolist()
        self.means = input_data['mean_psi']
        self.conf_interval = input_data['conf_interval']
        self.quartiles = input_data['quartiles']
        self.coords = input_data['coords']
        self.matrix_link = input_data['matrix_link']
        self.variances = input_data['variance_list']

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

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)