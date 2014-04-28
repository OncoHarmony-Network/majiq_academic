from __future__ import division


class Event(object):
    # TODO: Encapsulate calls regarding the mean, the conf_interval and the quartiles here

    def __init__(self, input_data=None):
        # Data for event table
        self.number = input_data['number']
        self.id = input_data['id']
        self.type = input_data['type']
        self.bins = input_data['bins']
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