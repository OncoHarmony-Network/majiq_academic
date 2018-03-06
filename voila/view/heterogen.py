import os
from tempfile import NamedTemporaryFile

from voila import constants
from voila.api import SpliceGraph
from voila.exceptions import NoLsvsFound, NotHeterogenVoilaFile
from voila.utils.voila_log import voila_log
from voila.view.html import Html


class Heterogen(Html):
    def __init__(self, args):
        super().__init__(args)
        with ViewHeterogen(args.voila_file) as m:
            if m.analysis_type != constants.ANALYSIS_PSI:
                raise NotHeterogenVoilaFile(args.voila_file)
            self.metadata = ViewHeterogenMatrix(m).metadata

        log = voila_log()
        # if not args.no_tsv:
        #     io_voila.het_tab_output(args)

        if not args.no_html:
            log.warning('No html output.')
        #     with Voila(args.voila_file, 'r') as v:
        #         self.metainfo = v.get_metainfo()
        #     self.render_summaries()
        #     self.render_index()
        #     copy_static(args)

    # @classmethod
    # def arg_parents(cls):
    #     return (
    #         cls.base_args(), cls.html_args(), cls.voila_file_args(), cls.multiproccess_args(), cls.output_args(),
    #         cls.lsv_id_search_args(), cls.gene_search_args()
    #     )

    def render_index(self):
        log = voila_log()
        log.info('Render Heterogen HTML index')
        log.debug('Start index render')
        args = self.args
        env = self.env
        metainfo = self.metainfo
        index_row_template = env.get_template('het_index_row.html')

        with NamedTemporaryFile(dir=self.get_template_dir()) as tmp_index_file:
            lsv_count = None

            # with Voila(args.voila_file, 'r') as v:
            #
            #     lsv_count = v.get_lsv_count(args)
            #     too_many_lsvs = lsv_count > constants.MAX_LSVS_HET_INDEX
            #
            #     for index, lsv in enumerate(v.get_voila_lsvs(args)):
            #         log.debug('Writing {0} to index'.format(lsv.lsv_id))
            #
            #         for el in index_row_template.generate(
            #                 lsv=lsv,
            #                 index=index,
            #                 link=self.voila_links[lsv.name],
            #                 too_many_lsvs=too_many_lsvs,
            #                 threshold=.2,
            #                 lexps=metainfo
            #         ):
            #             tmp_index_file.write(bytearray(el, encoding='utf-8'))

            if not tmp_index_file.tell():
                raise NoLsvsFound()

            log.debug('Write tmp index to actual index')

            with open(os.path.join(args.output, 'index.html'), 'w') as html:
                index_template = env.get_template('index_het_summary_template.html')
                for el in index_template.generate(
                        lexps=metainfo,
                        tmp_index_file_name=os.path.basename(tmp_index_file.name),
                        table_marks=self.table_marks_set(lsv_count),
                        lsvs_count=lsv_count,
                        prev_page=None,
                        next_page=None
                ):
                    html.write(el)

                log.debug('End index render')

    def render_summaries(self):
        log = voila_log()
        log.info('Render Heterogen HTML summaries')
        log.debug('Start summaries render')

        summary_template = self.env.get_template('het_summary_template.html')
        args = self.args
        summaries_subfolder = self.get_summaries_subfolder()
        metainfo = self.metainfo

        with SpliceGraph(args.splice_graph, 'r') as sg:
            gene_experiments_list = sg.get_experiments()
            prev_page = None
            page_count = sg.get_page_count(args)

            log.debug('There will be {0} pages of gene summaries'.format(page_count))

            for index, (lsv_dict, genes) in enumerate(sg.get_paginated_genes_with_lsvs(args)):

                # todo: fix this... this shouldn't be this way
                for gene_id, lsvs in lsv_dict.items():
                    for lsv in lsvs:
                        lsv.box_plot_data = [g.median.tolist() for g in lsv.het.groups]

                page_name = self.get_page_name(index)
                table_marks = tuple(self.table_marks_set(len(gene_set)) for gene_set in lsv_dict)
                next_page = self.get_next_page(index, page_count)

                experiments = tuple(
                    self.gene_experiments(exps, genes, gene_experiments_list) for exps in metainfo['experiment_names'])

                self.add_to_voila_links(lsv_dict, page_name)

                log.debug('Write page {0}'.format(page_name))

                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:

                    for el in summary_template.generate(
                            page_name=page_name,
                            genes_dict=lsv_dict,
                            genes_exps_list=experiments,
                            lexps=metainfo,
                            threshold=.2,
                            lsv_text_version=constants.LSV_TEXT_VERSION,
                            table_marks=table_marks,
                            prev_page=prev_page,
                            next_page=next_page,
                            gtf=args.gtf
                    ):
                        html.write(el)

                    prev_page = page_name

                    log.debug('End summaries render')
