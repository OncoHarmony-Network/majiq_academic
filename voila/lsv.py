from __future__ import division


class Lsv(object):

    def __init__(self, input_data=None):

        self.id = input_data['id']
        self.number = input_data['number']
        self.type = input_data['type']
        self.bins = input_data['bins_list']
        self.mean_psi = input_data['mean_psi']
        self.conf_interval = input_data['conf_interval']
        self.quartiles = input_data['quartiles']
        self.coords = input_data['coords']
        self.matrix_link = input_data['matrix_link']

    def get_bins(self):
        return self.bins

    def get_id(self):
        return self.id

    def set_excl_incl(self, excl_incl_set):
        self.excl_incl = excl_incl_set

    def get_excl_incl(self):
        return self.excl_incl

    def get_matrix_link(self):
        return self.matrix_link