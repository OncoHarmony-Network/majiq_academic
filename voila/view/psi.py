import errno
import os

from voila import constants, io_voila
from voila.api.view_matrix import ViewMatrix
from voila.api.view_splice_graph import PsiSpliceGraph
from voila.utils.run_voila_utils import table_marks_set, copy_static, get_env
from voila.utils.voila_log import voila_log
from voila.utils.voila_pool import VoilaPool
from voila.view.html import Html
from voila.voila_args import VoilaArgs


def create_gene_db(gene_ids, args, experiment_names):
    env = get_env()
    log = voila_log()
    with PsiSpliceGraph(args.splice_graph) as sg, ViewMatrix(args.voila_file) as m:
        for gene_id in gene_ids:
            log.debug('creating {}'.format(gene_id))
            with open(os.path.join(args.output, 'db', '{}.js'.format(gene_id)), 'w') as f:
                for el in env.get_template('gene_db_template.html').generate(
                        gene=sg.gene(gene_id).get.get_experiment(experiment_names),
                        lsvs=(m.psi(lsv_id) for lsv_id in m.lsv_ids(gene_id))
                ):
                    f.write(el)


class Psi(Html, VoilaArgs):
    def __init__(self, args):
        """
        Render psi output.
        :param args: command line arguments
        :return: None
        """
        super(Psi, self).__init__(args)

        # if not args.no_html:
        #     with ViewMatrix(args.voila_file, 'r') as m:
        #         self.metadatda = m.metadata
        #     self.render_summaries()
        #     self.render_index()
        #     copy_static(args)

        if not args.no_tsv:
            io_voila.tab_output(args, self.voila_links)

        # if args.gtf:
        #     io_voila.generic_feature_format_txt_files(args)
        #
        # if args.gff:
        #     io_voila.generic_feature_format_txt_files(args, out_gff3=True)

    def render_index(self):
        log = voila_log()
        log.info('Render PSI HTML index')
        log.debug('Start index render')

        args = self.args
        metadata = self.metadatda

        with ViewMatrix(args.voila_file, 'r') as m:
            lsv_count = m.get_lsv_count(args)
            lsv_ids = tuple(m.get_lsvs(args))
            gene_ids = tuple(m.gene_ids)

        too_many_lsvs = lsv_count > constants.MAX_LSVS_PSI_INDEX

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
            for lsvs, genes in zip(chunkify(lsv_ids, vp.processes), chunkify(gene_ids, vp.processes)):
                multiple_results.append(vp.pool.apply_async(create_gene_db, (genes, args, names)))

            for res in multiple_results:
                res.get()

        log.debug('finished writing db files.')

        with ViewMatrix(args.voila_file, 'r') as m, PsiSpliceGraph(args.splice_graph) as sg:
            with open(os.path.join(args.output, 'index.html'), 'w') as html:
                index_template = self.env.get_template('index_single_summary_template.html')
                for el in index_template.generate(
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
                        lsv_text_version=constants.LSV_TEXT_VERSION,
                ):
                    html.write(el)

    @classmethod
    def arg_parents(cls):
        return (cls.base_args(), cls.output_args(), cls.html_args(), cls.gene_search_args(), cls.lsv_type_search_args(),
                cls.lsv_id_search_args(), cls.voila_file_args())

    def render_summaries(self):
        log = voila_log()
        log.info('Render PSI HTML summaries')
        log.debug('Start summaries render')

        summary_template = self.env.get_template("psi_summary_template.html")
        args = self.args
        summaries_subfolder = self.get_summaries_subfolder()
        metadata = self.metadatda
        database_name = self.database_name()

        # with Pool(processes=4) as pool:

        with PsiSpliceGraph(args.splice_graph, 'r') as sg:
            genome = sg.genome
            prev_page = None
            page_count = sg.get_page_count(args)

            for index, (lsv_dict, genes) in enumerate(sg.get_paginated_genes_with_lsvs(args)):
                page_name = self.get_page_name(index)
                table_marks = tuple(table_marks_set(len(gene_set)) for gene_set in lsv_dict)
                next_page = self.get_next_page(index, page_count)

                self.add_to_voila_links(lsv_dict, page_name)

                log.debug('Writing {0}'.format(page_name))
                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:

                    for el in summary_template.generate(
                            genes=[sg.gene(gene_id) for gene_id in genes],
                            table_marks=table_marks,
                            lsvs=lsv_dict,
                            prev_page=prev_page,
                            next_page=next_page,
                            namePage=page_name,
                            metadata=metadata,
                            lsv_text_version=constants.LSV_TEXT_VERSION,
                            gtf=args.gtf,
                            database_name=database_name,
                            genome=genome
                    ):
                        html.write(el)

                prev_page = page_name

        log.debug('End summaries render')
