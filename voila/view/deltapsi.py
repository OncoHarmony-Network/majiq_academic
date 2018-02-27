import errno
import os

from voila import constants, io_voila
from voila.api.matrix_hdf5 import lsv_id_to_gene_id
from voila.api.view_matrix import ViewDeltaPsi, ViewDeltaPsiMatrix
from voila.api.view_splice_graph import ViewSpliceGraph, ViewGene
from voila.exceptions import NotDeltaPsiVoilaFile
from voila.utils.voila_log import voila_log
from voila.utils.voila_pool import VoilaPool
from voila.view.html import Html


class DeltaPsi(Html):
    def __init__(self, args):
        super(DeltaPsi, self).__init__(args)

        if not args.disable_html:
            with ViewDeltaPsi(args.voila_file) as m:
                if m.analysis_type != constants.ANALYSIS_DELTAPSI:
                    raise NotDeltaPsiVoilaFile(args.voila_file)
                self.metadata = m.metadata

            self.copy_static(args)
            self.create_db_files()
            self.render_summaries()
            self.render_index()

        if not args.disable_tsv:
            io_voila.delta_psi_tab_output(args, self.voila_links)

    def render_index(self):
        log = voila_log()
        log.info('Render Delta PSI HTML index')
        log.debug('Start index render')
        args = self.args
        metadata = self.metadata

        with ViewSpliceGraph(args.splice_graph) as sg, ViewDeltaPsi(args.voila_file) as m:
            lsv_count = ViewDeltaPsiMatrix(m).view_lsv_count(args)
            too_many_lsvs = lsv_count > constants.MAX_LSVS_DELTAPSI_INDEX

            with open(os.path.join(args.output, 'index.html'), 'w') as html:
                index_template = self.get_env().get_template('index_delta_summary_template.html')
                html.write(
                    index_template.render(
                        lexps=metadata,
                        table_marks=self.table_marks_set(lsv_count),
                        lsvs_count=lsv_count,
                        prev_page=None,
                        next_page=None,
                        database_name='deltapsi_' + self.database_name(),
                        lsvs=list(m.delta_psi(lsv_id) for lsv_id in ViewDeltaPsiMatrix(m).view_lsv_ids(args)),
                        genes=sg.gene,
                        gene_ids=set(lsv_id_to_gene_id(lsv_id) for lsv_id in ViewDeltaPsiMatrix(m).view_lsv_ids(args)),
                        links=self.voila_links,
                        too_many_lsvs=too_many_lsvs,
                        metadata=metadata,
                        genome=sg.genome,
                        lsv_text_version=constants.LSV_TEXT_VERSION,
                        threshold=args.threshold
                    )
                )

                log.debug('End index render')

    @classmethod
    def create_gene_db(cls, gene_ids, args, experiment_names):
        template = cls.get_env().get_template('gene_db_template.html')
        log = voila_log()
        with ViewSpliceGraph(args.splice_graph) as sg, ViewDeltaPsi(args.voila_file) as m:
            for gene_id in gene_ids:
                log.debug('creating {}'.format(gene_id))

                with open(os.path.join(args.output, 'db', '{}.js'.format(gene_id)), 'w') as html:
                    html.write(
                        template.render(
                            gene=ViewGene(sg.gene(gene_id).get).get_experiment(experiment_names),
                            lsvs=tuple(
                                m.delta_psi(lsv_id) for lsv_id in ViewDeltaPsiMatrix(m).view_gene_lsvs(args, gene_id))
                        )
                    )

    @classmethod
    def create_summary(cls, metadata, args, database_name, paged):
        summary_template = cls.get_env().get_template("deltapsi_summary_template.html")
        summaries_subfolder = cls.get_summaries_subfolder(args)
        group_names = metadata['group_names']
        links = {}

        with ViewSpliceGraph(args.splice_graph, 'r') as sg, ViewDeltaPsi(args.voila_file) as m:
            genome = sg.genome
            page_count = ViewDeltaPsiMatrix(m).page_count(args)

            for index, genes in paged:
                page_name = cls.get_page_name(args, index)
                next_page = cls.get_next_page(args, index, page_count)
                prev_page = cls.get_prev_page(args, index)

                lsv_dict = {gene_id: tuple(lsv_id for lsv_id in ViewDeltaPsiMatrix(m).view_gene_lsvs(args, gene_id)) for
                            gene_id in genes}
                table_marks = tuple(cls.table_marks_set(len(gene_set)) for gene_set in lsv_dict)

                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:
                    html.write(
                        summary_template.render(
                            page_name=page_name,
                            threshold=args.threshold,
                            lsv_text_version=constants.LSV_TEXT_VERSION,
                            table_marks=table_marks,
                            prev_page=prev_page,
                            next_page=next_page,
                            gtf=args.gtf,
                            group_names=group_names,
                            genes=list(sg.gene(gene_id).get for gene_id in genes),
                            lsv_ids=lsv_dict,
                            delta_psi_lsv=m.delta_psi,
                            metadata=metadata,
                            database_name=database_name,
                            genome=genome
                        )
                    )

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
            paged_genes = tuple(ViewDeltaPsiMatrix(m).paginated_genes(args))

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
            gene_ids = tuple(ViewDeltaPsiMatrix(m).view_gene_ids(args))

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
