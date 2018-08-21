import csv
import os
import uuid
from multiprocessing import Lock

from voila import constants
from voila.api.view_matrix import ViewPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.exceptions import NotPsiVoilaFile
from voila.utils.voila_log import voila_log
from voila.view.html import Html
from voila.view.tsv import Tsv

lock = Lock()


class Psi(Html, Tsv):
    def __init__(self, args):
        super().__init__(args, ViewPsi)
        self.db_id = uuid.uuid4().hex

        with ViewPsi(args) as m:
            if m.analysis_type != constants.ANALYSIS_PSI:
                raise NotPsiVoilaFile(args)
            self.view_metadata = m.view_metadata

        if not args.disable_html:
            self.copy_static()
            if not args.disable_db:
                self.render_dbs()
            self.render_html('psi_index.html', 'psi_summary.html')

        if not args.disable_tsv:
            self.psi_tab_output()

    def tsv_row(self, gene_ids, tsv_file, fieldnames):
        args = self.args
        voila_links = self.voila_links
        log = voila_log()
        with ViewPsi(args) as m, ViewSpliceGraph(args) as sg:

            with open(tsv_file, 'a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene_id in gene_ids:
                    gene = sg.gene(gene_id)
                    for lsv_id in m.view_gene_lsvs(gene.id):
                        lsv = m.lsv(lsv_id)
                        lsv_junctions = lsv.junctions
                        annot_juncs = sg.annotated_junctions(gene, lsv)
                        lsv_exons = sg.lsv_exons(gene, lsv)

                        row = {
                            '#Gene Name': gene.name,
                            'Gene ID': gene.id,
                            'LSV ID': lsv_id,
                            'LSV Type': lsv.lsv_type,
                            'A5SS': lsv.a5ss,
                            'A3SS': lsv.a3ss,
                            'ES': lsv.exon_skipping,
                            'Num. Junctions': lsv.junction_count,
                            'Num. Exons': lsv.exon_count,
                            'chr': gene.chromosome,
                            'strand': gene.strand,
                            'De Novo Junctions': self.semicolon_join(annot_juncs),
                            'Junctions coords': self.semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'Exons coords': self.semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in self.filter_exons(lsv_exons)
                            ),
                            # 'IR coords': self.semicolon_join(
                            #     '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                            # ),
                            'E(PSI) per LSV junction': self.semicolon_join(lsv.means),
                            'Var(E(PSI)) per LSV junction': self.semicolon_join(lsv.variances)
                        }

                        if voila_links:
                            summary_path = voila_links[gene.id]
                            if not os.path.isabs(summary_path):
                                summary_path = os.path.join(os.getcwd(), args.output, summary_path)
                            row['Voila link'] = "file://{0}".format(summary_path)

                        log.debug('Write TSV row for {0}'.format(lsv_id))

                        lock.acquire()
                        writer.writerow(row)
                        lock.release()

    def psi_tab_output(self):
        voila_links = self.voila_links

        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction',
                      'LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords']
        if voila_links:
            fieldnames.append('Voila link')

        self.write_tsv(fieldnames, ViewPsi)
