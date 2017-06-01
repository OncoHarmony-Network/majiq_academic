from voila import io_voila
from voila.io_voila import Voila
from voila.voila_args import VoilaArgs


class Heterogen(VoilaArgs):
    def __init__(self, args):
        with Voila(args.voila_file, 'r') as v:
            metainfo = v.get_metainfo()
            lsvs = v.get_voila_lsvs()

        if not args.no_tsv:
            io_voila.het_tab_output(args, lsvs, metainfo)

    @classmethod
    def arg_parents(cls):
        return (
            cls.base_args(), cls.html_args(), cls.voila_file_args(), cls.multiproccess_args(), cls.output_args()
        )
