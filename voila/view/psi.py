import errno
import os

from voila import constants, io_voila
from voila.api.view_matrix import ViewPsi
from voila.api.view_splice_graph import ViewSpliceGraph, ViewGene
from voila.utils.exceptions import NotPsiVoilaFile
from voila.utils.run_voila_utils import table_marks_set, copy_static, get_env
from voila.utils.voila_log import voila_log
from voila.utils.voila_pool import VoilaPool
from voila.view.html import Html


class Psi(Html):
    def __init__(self, args):
        """
        Render psi output.
        :param args: command line arguments
        :return: None
        """
        super(Psi, self).__init__(args)

        if not args.disable_html:
            with ViewPsi(args.voila_file) as m:
                if m.analysis_type != constants.ANALYSIS_PSI:
                    raise NotPsiVoilaFile(args.voila_file)
                self.metadata = m.metadata
            copy_static(args)
            self.create_db_files()
            self.render_summaries()
            self.render_index()

        if not args.disable_tsv:
            io_voila.psi_tab_output(args, self.voila_links)

    def create_db_files(self):
        args = self.args
        metadata = self.metadata
        log = voila_log()
        log.info('Create DB files')
        with ViewPsi(args.voila_file, 'r') as m:
            gene_ids = tuple(m.gene_ids)

        try:
            os.makedirs(os.path.join(args.output, 'db'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        names = metadata['experiment_names']
        multiple_results = []

        with VoilaPool() as vp:
            for genes in self.chunkify(gene_ids, vp.processes):
                multiple_results.append(vp.pool.apply_async(self.create_gene_db, (genes, args, names)))

            for res in multiple_results:
                res.get()

        log.debug('finished writing db files.')

    def render_index(self):
        log = voila_log()
        log.info('Render PSI HTML index')
        log.debug('Start index render')

        args = self.args
        metadata = self.metadata

        with ViewPsi(args.voila_file) as m, ViewSpliceGraph(args.splice_graph) as sg:
            lsv_ids = tuple(m.view_lsv_ids(args))
            lsv_count = m.view_lsv_count(args)
            too_many_lsvs = lsv_count > constants.MAX_LSVS_PSI_INDEX

            with open(os.path.join(args.output, 'index.html'), 'w') as html:
                index_template = self.env.get_template('index_single_summary_template.html')
                html.write(
                    index_template.render(
                        lexps=metadata,
                        metadata=metadata,
                        lsvs_count=lsv_count,
                        table_marks=table_marks_set(lsv_count),
                        prev_page=None,
                        next_page=None,
                        links=self.voila_links,
                        genes=sg.gene,
                        lsvs=(m.psi(lsv_id) for lsv_id in lsv_ids),
                        too_many_lsvs=too_many_lsvs,
                        database_name=self.database_name(),
                        genome=sg.genome,
                        lsv_text_version=constants.LSV_TEXT_VERSION
                    )
                )

    @classmethod
    def create_summary(cls, metadata, args, database_name, paged):
        summary_template = get_env().get_template("psi_summary_template.html")
        summaries_subfolder = cls.get_summaries_subfolder(args)
        links = {}

        with ViewSpliceGraph(args.splice_graph, 'r') as sg, ViewPsi(args.voila_file) as m:
            genome = sg.genome
            page_count = m.page_count(args)

            for index, genes in paged:
                page_name = cls.get_page_name(args, index)
                next_page = cls.get_next_page(args, index, page_count)
                prev_page = cls.get_prev_page(args, index)
                lsv_dict = {gene_id: tuple(lsv_id for lsv_id in m.view_lsv_ids(args, gene_id)) for gene_id in genes}
                table_marks = tuple(table_marks_set(len(gene_set)) for gene_set in lsv_dict)

                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:
                    html.write(
                        summary_template.render(
                            genes=list(sg.gene(gene_id).get for gene_id in genes),
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

                links.update(dict(cls.voila_links(lsv_dict, page_name)))
            return links

    def render_summaries(self):
        log = voila_log()
        log.info('Render PSI HTML summaries')
        log.debug('Start summaries render')

        args = self.args
        metadata = self.metadata
        database_name = self.database_name()

        with ViewPsi(args.voila_file) as m:
            paged_genes = tuple(m.paginated_genes(args))

        multiple_results = []
        with VoilaPool() as vp:
            for paged in self.chunkify(tuple(enumerate(paged_genes)), vp.processes):
                multiple_results.append(
                    vp.pool.apply_async(self.create_summary, (metadata, args, database_name, paged)))

            for res in multiple_results:
                self.voila_links.update(res.get())

        log.debug('End summaries render')

    @staticmethod
    def create_gene_db(gene_ids, args, experiment_names):
        template = get_env().get_template('gene_db_template.html')
        with ViewSpliceGraph(args.splice_graph) as sg, ViewPsi(args.voila_file) as m:
            for gene_id in gene_ids:
                with open(os.path.join(args.output, 'db', '{}.js'.format(gene_id)), 'w') as html:
                    html.write(
                        template.render(
                            gene=ViewGene(sg.gene(gene_id).get).get_experiment(experiment_names),
                            lsvs=(m.psi(lsv_id) for lsv_id in m.view_lsv_ids(args, gene_id))
                        )
                    )
