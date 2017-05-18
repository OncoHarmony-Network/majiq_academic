from voila import io_voila
from voila.io_voila import Voila


class Heterogen(object):
    def __init__(self, args):
        with Voila(args.voila_file, 'r') as v:
            metainfo = v.get_metainfo()
            lsvs = v.get_voila_lsvs()

        if not args.no_tsv:
            io_voila.het_tab_output(args, lsvs, metainfo)
