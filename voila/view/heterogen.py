from voila import io_voila
from voila.voila_args import VoilaArgs


class Heterogen(VoilaArgs):
    def __init__(self, args):
        if not args.no_tsv:
            io_voila.het_tab_output(args)

    @classmethod
    def arg_parents(cls):
        return (
            cls.base_args(), cls.html_args(), cls.voila_file_args(), cls.multiproccess_args(), cls.output_args(), cls.lsv_id_search_args(), cls.gene_search_args()
        )
