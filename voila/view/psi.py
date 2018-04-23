import csv
import os
from multiprocessing import Lock

from voila import constants
from voila.api.matrix_hdf5 import lsv_id_to_gene_id
from voila.api.view_matrix import ViewPsi
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.exceptions import NotPsiVoilaFile
from voila.utils.voila_log import voila_log
from voila.view.html import Html
from voila.view.tsv import Tsv

lock = Lock()


class Psi(Html, Tsv):
    def __init__(self, args):
        """
        Render psi output.
        :param args: command line arguments
        :return: None
        """
        super(Psi, self).__init__(args)

        with ViewPsi(args) as m:
            if m.analysis_type != constants.ANALYSIS_PSI:
                raise NotPsiVoilaFile(args.voila_file)
            self.view_metadata = m.view_metadata

        if not args.disable_html:
            self.copy_static()
            self.render_dbs()
            self.render_summaries()
            self.render_index()

        if not args.disable_tsv:
            self.psi_tab_output()

    def render_index(self):
        log = voila_log()
        log.info('Render PSI HTML index')
        log.debug('Start index render')

        args = self.args
        metadata = self.view_metadata

        with ViewPsi(args) as m, ViewSpliceGraph(args) as sg:
            lsv_ids = tuple(m.view_lsv_ids())
            lsv_count = m.view_lsv_count()
            too_many_lsvs = lsv_count > constants.MAX_LSVS_PSI_INDEX

            with open(os.path.join(args.output, 'index.html'), 'w') as html:
                index_template = self.get_env().get_template('index_single_summary_template.html')
                html.write(
                    index_template.render(
                        lexps=metadata,
                        metadata=metadata,
                        lsvs_count=lsv_count,
                        table_marks=self.table_marks_set(lsv_count),
                        prev_page=None,
                        next_page=None,
                        links=self.voila_links,
                        genes=sg.gene,
                        gene_ids=set(lsv_id_to_gene_id(lsv_id) for lsv_id in m.view_lsv_ids()),
                        lsvs=(m.psi(lsv_id) for lsv_id in lsv_ids),
                        too_many_lsvs=too_many_lsvs,
                        database_name=self.database_name(),
                        genome=sg.genome,
                        lsv_text_version=constants.LSV_TEXT_VERSION
                    )
                )

    def create_summary(self, paged):
        args = self.args
        summary_template = self.get_env().get_template("psi_summary_template.html")
        summaries_subfolder = self.get_summaries_subfolder(args)
        database_name = self.database_name()
        metadata = self.view_metadata
        links = {}

        with ViewSpliceGraph(args) as sg, ViewPsi(args) as m:
            genome = sg.genome
            page_count = m.page_count()

            for index, genes in paged:
                page_name = self.get_page_name(args, index)
                next_page = self.get_next_page(args, index, page_count)
                prev_page = self.get_prev_page(args, index)
                lsv_dict = {gene_id: tuple(lsv_id for lsv_id in m.view_gene_lsvs(gene_id)) for
                            gene_id in genes}
                table_marks = tuple(self.table_marks_set(len(gene_set)) for gene_set in lsv_dict)

                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:
                    html.write(
                        summary_template.render(
                            genes=list(sg.gene(gene_id) for gene_id in genes),
                            table_marks=table_marks,
                            lsv_ids=lsv_dict,
                            psi_lsv=m.psi,
                            prev_page=prev_page,
                            next_page=next_page,
                            namePage=page_name,
                            metadata=metadata,
                            lsv_text_version=constants.LSV_TEXT_VERSION,
                            gtf=args.gtf,
                            database_name=database_name,
                            genome=genome
                        )
                    )

                links.update(self.get_voila_links(lsv_dict, page_name))

            return links

    def render_summaries(self):
        # log = voila_log()
        # log.info('Render PSI HTML summaries')
        # log.debug('Start summaries render')
        #
        # args = self.args
        # metadata = self.view_metadata
        # database_name = self.database_name()
        #
        # with ViewPsi(args) as m:
        #     paged_genes = tuple(m.paginated_genes())
        #
        # multiple_results = []
        # with VoilaPool() as vp:
        #     for paged in self.chunkify(list(enumerate(paged_genes)), vp.processes):
        #         multiple_results.append(
        #             vp.apply_async(self.create_summary, (metadata, args, database_name, paged)))
        #
        #     for res in multiple_results:
        #         self.voila_links.update(res.get())
        #
        # log.debug('End summaries render')
        self.create_summaries(ViewPsi)

    def render_dbs(self):
        self.create_db_files(ViewPsi, 'psi')

    def tsv_row(self, gene_ids, tsv_file, fieldnames):
        args = self.args
        voila_links = self.voila_links
        log = voila_log()
        with ViewPsi(args) as m, ViewSpliceGraph(args) as sg:

            with open(tsv_file, 'a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene in sg.genes(gene_ids):
                    for lsv_id in m.view_gene_lsvs(gene.id):
                        lsv = m.psi(lsv_id)
                        lsv_junctions = list(gene.lsv_junctions(lsv))
                        lsv_exons = list(gene.lsv_exons(lsv))

                        row = {
                            '#Gene Name': gene.name,
                            'Gene ID': gene.id,
                            'LSV ID': lsv_id,
                            'LSV Type': lsv.lsv_type,
                            'A5SS': lsv.prime5,
                            'A3SS': lsv.prime3,
                            'ES': lsv.exon_skipping,
                            'Num. Junctions': lsv.junction_count,
                            'Num. Exons': lsv.exon_count,
                            'chr': gene.chromosome,
                            'strand': gene.strand,
                            'De Novo Junctions': self.semicolon_join(
                                int(not junc.annotated) for junc in lsv_junctions
                            ),
                            'Junctions coords': self.semicolon_join(
                                '{0}-{1}'.format(junc.start, junc.end) for junc in lsv_junctions
                            ),
                            'Exons coords': self.semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in self.filter_exons(lsv_exons)
                            ),
                            'IR coords': self.semicolon_join(
                                '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                            ),
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

        self.write_tsv(fieldnames)
