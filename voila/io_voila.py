__author__ = 'abarrera'
import cPickle as pkl


class VoilaInput(object):
    """Standard input interface by experiment used by Voila"""
    def __init__(self, lsvs=(), metainfo=None):
        self.lsvs = lsvs
        self.metainfo = metainfo

    def get_lsvs(self):
        return self.lsvs

    def add_lsv(self, l):
        self.lsvs.append(l)

    def get_metainfo(self):
        """Retrieve experiment metainfo (primarily samples information, expandable)"""
        return self.metainfo

    def samples_metainfo(self):
        """Sample information generator, one by one."""
        for sample_info in self.metainfo:
            yield sample_info


def dump_voila_input(voila_input, target, logger=None, protocol=-1):
    import os
    base_path, name = os.path.split(target)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    try:
        pkl.dump(voila_input, target, protocol)
    except pkl.PickleError, e:
        if logger:
            logger.error("Dumping file %s in %s:\n\t%s." % (voila_input, target, e.message), exc_info=1)


def load_voila_input(voila_input_file, logger=None):
    try:
        with open(voila_input_file, 'rb') as vif:
            voila_input = pkl.load(vif)
            return voila_input
    except pkl.PickleError, e:
        if logger:
            logger.error("Loading the file %s:\n\t%s." % (voila_input_file, e.message), exc_info=1)