import csv
import os
from abc import ABC

from voila.api.view_matrix import ViewHeterogens
from voila.utils.voila_log import voila_log
from voila.processes import VoilaPool, VoilaQueue
from voila.view.html import Html


class Tsv(ABC):
    def __init__(self, args):
        self.args = args

    @staticmethod
    def semicolon_join(value_list):
        return ';'.join(str(x) for x in value_list)

    @staticmethod
    def filter_exons(exons):
        for exon in exons:
            if exon.start == -1:
                yield 'nan', exon.end
            elif exon.end == -1:
                yield exon.start, 'nan'
            else:
                yield exon.start, exon.end

    def write_tsv(self, fieldnames):
        log = voila_log()
        log.info("Creating Tab-delimited output file")

        args = self.args
        output_html = Html.get_output_html(args, args.voila_file[0])
        tsv_file = os.path.join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

        with ViewHeterogens(args) as m:
            view_gene_ids = list(m.gene_ids)

        with open(tsv_file, 'w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

        with VoilaPool() as vp, VoilaQueue() as vq:
            for gene_ids in self.chunkify(view_gene_ids, vp.processes):
                multiple_results.append(vp.apply_async(self.tsv_row, (gene_ids, tsv_file, fieldnames)))

        [r.get() for r in multiple_results]

        log.info("Delimited output file successfully created in: %s" % tsv_file)
