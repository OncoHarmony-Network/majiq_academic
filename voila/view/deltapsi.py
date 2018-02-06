import errno
import os

from voila import constants, io_voila
from voila.api.view_matrix import ViewDeltaPsi
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.utils.run_voila_utils import table_marks_set, copy_static, get_env
from voila.utils.voila_log import voila_log
from voila.utils.voila_pool import VoilaPool
from voila.view.html import Html


class DeltaPsi(Html):
    def __init__(self, args):
        super(DeltaPsi, self).__init__(args)

        if not args.disable_html:
            copy_static(args)
            with ViewDeltaPsi(args.voila_file) as m:
                self.metadata = m.metadata
                self.lsv_ids = tuple(m.view_lsv_ids(args))
            self.create_db_files()
            self.render_summaries()
            self.render_index()

        if not args.disable_tsv:
            io_voila.delta_psi_tab_output(args, self.voila_links)

        # if args.gtf:
        #     io_voila.generic_feature_format_txt_files(args)
        #
        # if args.gff:
        #     io_voila.generic_feature_format_txt_files(args, out_gff3=True)

    def render_index(self):
        log = voila_log()
        log.info('Render Delta PSI HTML index')
        log.debug('Start index render')
        args = self.args
        env = self.env
        metadata = self.metadata

        with ViewSpliceGraph(args.splice_graph) as sg, ViewDeltaPsi(args.voila_file) as m:
            lsv_count = m.view_lsv_count(args)
            too_many_lsvs = lsv_count > constants.MAX_LSVS_DELTAPSI_INDEX

            with open(os.path.join(args.output, 'index.html'), 'w') as html:
                index_template = env.get_template('index_delta_summary_template.html')
                for el in index_template.generate(
                        lexps=metadata,
                        table_marks=table_marks_set(lsv_count),
                        lsvs_count=lsv_count,
                        prev_page=None,
                        next_page=None,
                        database_name=self.database_name(),
                        lsvs=(m.delta_psi(lsv_id) for lsv_id in self.lsv_ids),
                        genes=sg.gene,
                        links=self.voila_links,
                        too_many_lsvs=too_many_lsvs,
                        metadata=metadata,
                        genome=sg.genome,
                        lsv_text_version=constants.LSV_TEXT_VERSION,
                        threshold=args.threshold
                ):
                    html.write(el)

                log.debug('End index render')

    @staticmethod
    def create_gene_db(gene_ids, args, experiment_names):
        env = get_env()
        log = voila_log()
        with ViewSpliceGraph(args.splice_graph) as sg, ViewDeltaPsi(args.voila_file) as m:
            for gene_id in gene_ids:
                log.debug('creating {}'.format(gene_id))
                with open(os.path.join(args.output, 'db', '{}.js'.format(gene_id)), 'w') as f:
                    for el in env.get_template('gene_db_template.html').generate(
                            gene=sg.gene(gene_id).get.get_experiment(experiment_names),
                            lsvs=(m.delta_psi(lsv_id) for lsv_id in m.lsv_ids(gene_id))
                    ):
                        f.write(el)

    @classmethod
    def create_summary(cls, metadata, args, database_name, paged):

        summary_template = get_env().get_template("deltapsi_summary_template.html")
        summaries_subfolder = cls.get_summaries_subfolder(args)
        group_names = metadata['group_names']
        links = {}

        with ViewSpliceGraph(args.splice_graph, 'r') as sg, ViewDeltaPsi(args.voila_file) as m:
            genome = sg.genome
            page_count = m.get_page_count(args)

            for index, genes in paged:
                page_name = cls.get_page_name(args, index)
                next_page = cls.get_next_page(args, index, page_count)
                prev_page = cls.get_prev_page(args, index)
                lsv_dict = {gene_id: tuple(lsv_id for lsv_id in m.view_lsv_ids(args, gene_id)) for gene_id in genes}
                table_marks = tuple(table_marks_set(len(gene_set)) for gene_set in lsv_dict)

                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:
                    for el in summary_template.generate(
                            page_name=page_name,
                            threshold=args.threshold,
                            lsv_text_version=constants.LSV_TEXT_VERSION,
                            table_marks=table_marks,
                            prev_page=prev_page,
                            next_page=next_page,
                            gtf=args.gtf,
                            group_names=group_names,
                            genes=[sg.gene(gene_id) for gene_id in genes],
                            lsv_ids=lsv_dict,
                            delta_psi_lsv=m.delta_psi,
                            metadata=metadata,
                            database_name=database_name,
                            genome=genome
                    ):
                        html.write(el)

                links.update(dict(cls.voila_links(lsv_dict, page_name)))

            return links

    def render_summaries(self):
        log = voila_log()
        log.info('Render Delta PSI HTML summaries')
        log.debug('Start summaries render')

        args = self.args
        metadata = self.metadata
        database_name = self.database_name()

        with ViewDeltaPsi(args.voila_file) as m:
            paged_genes = tuple(m.paginated_genes(args))

        multiple_results = []
        with VoilaPool() as vp:
            for paged in self.chunkify(tuple(enumerate(paged_genes)), vp.processes):
                multiple_results.append(
                    vp.pool.apply_async(self.create_summary, (metadata, args, database_name, paged)))

            for res in multiple_results:
                self.voila_links.update(res.get())

    def create_db_files(self):
        args = self.args
        metadata = self.metadata
        log = voila_log()
        log.info('Create DB files')

        with ViewDeltaPsi(args.voila_file) as m:
            gene_ids = tuple(m.view_gene_ids(args))

        try:
            os.makedirs(os.path.join(args.output, 'db'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        def chunkify(lst, n):
            for i in range(n):
                yield lst[i::n]

        names = metadata['experiment_names']
        multiple_results = []

        with VoilaPool() as vp:
            for genes in chunkify(gene_ids, vp.processes):
                multiple_results.append(vp.pool.apply_async(self.create_gene_db, (genes, args, names)))

            for res in multiple_results:
                res.get()

        log.debug('finished writing db files.')
