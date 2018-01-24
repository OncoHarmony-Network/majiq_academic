import os
import tempfile

from voila import io_voila, constants
from voila.api import Voila
from voila.api.view_matrix import ViewMatrix
from voila.api.view_splice_graph import DeltaPsiSpliceGraph
from voila.utils.exceptions import NoLsvsFound
from voila.utils.run_voila_utils import table_marks_set, copy_static
from voila.utils.voila_log import voila_log
from voila.view.html import Html
from voila.voila_args import VoilaArgs


class Deltapsi(Html, VoilaArgs):
    def __init__(self, args):
        super(Deltapsi, self).__init__(args)

        if not args.no_html:
            copy_static(args)
            with ViewMatrix(args.voila_file) as m:
                self.metadata = m.metadata
            self.render_summaries()
            # self.render_index()

        # if not args.no_tsv:
        #     io_voila.tab_output(args, self.voila_links)

        if args.gtf:
            io_voila.generic_feature_format_txt_files(args)

        if args.gff:
            io_voila.generic_feature_format_txt_files(args, out_gff3=True)

    @classmethod
    def arg_parents(cls):
        # base, html, gene_search, lsv_type_search, lsv_id_search, voila_file,
        #                                parser_delta, multiprocess, output

        parser = cls.get_parser()
        # Probability threshold used to sum the accumulative probability of inclusion/exclusion.
        parser.add_argument('--threshold',
                            type=float,
                            default=0.2,
                            help='Filter out LSVs with no junction predicted to change over a certain value (in '
                                 'percentage).')

        parser.add_argument('--show-all',
                            dest='show_all',
                            action='store_true',
                            default=False,
                            help='Show all LSVs including those with no junction with significant change predicted.')

        return (
            cls.base_args(), cls.html_args(), cls.gene_search_args(), cls.lsv_type_search_args(),
            cls.lsv_id_search_args(), cls.voila_file_args(), cls.multiproccess_args(), cls.output_args(), parser
        )

    def render_index(self):
        log = voila_log()
        log.info('Render Delta PSI HTML index')
        log.debug('Start index render')
        args = self.args
        env = self.env
        metainfo = self.metainfo
        index_row_template = env.get_template('deltapsi_index_row.html')

        tmp_index_file = tempfile.mkstemp()[1]

        with open(tmp_index_file, 'w') as f:

            with Voila(args.voila_file) as v:

                lsv_count = v.get_lsv_count(args)
                too_many_lsvs = lsv_count > constants.MAX_LSVS_DELTAPSI_INDEX

                for index, lsv in enumerate(v.get_voila_lsvs(args)):
                    log.debug('Writing {0} to index'.format(lsv.lsv_id))

                    for el in index_row_template.generate(
                            lsv=lsv,
                            index=index,
                            link=self.voila_links[lsv.name],
                            too_many_lsvs=too_many_lsvs,
                            threshold=args.threshold,
                            lexps=metainfo
                    ):
                        # tmp_index_file.write(bytearray(el, encoding='utf-8'))
                        f.write(el)

            if not f.tell():
                raise NoLsvsFound()

        with open(tmp_index_file) as f:

            log.debug('Write tmp index to actual index')

            with open(os.path.join(args.output, 'index.html'), 'w') as html:
                index_template = env.get_template('index_delta_summary_template.html')
                for el in index_template.generate(
                        lexps=metainfo,
                        tmp_index_file=f.read(),
                        table_marks=table_marks_set(lsv_count),
                        lsvs_count=lsv_count,
                        prev_page=None,
                        next_page=None
                ):
                    html.write(el)

                log.debug('End index render')

        os.remove(tmp_index_file)

    def render_summaries(self):
        log = voila_log()
        log.info('Render Delta PSI HTML summaries')
        log.debug('Start summaries render')

        summary_template = self.env.get_template('deltapsi_summary_template.html')
        args = self.args
        summaries_subfolder = self.get_summaries_subfolder()
        metadata = self.metadata

        with DeltaPsiSpliceGraph(args.splice_graph) as sg:
            prev_page = None
            page_count = sg.get_page_count()
            genome = sg.genome
            database_name = self.database_name()
            log.debug('There will be {0} pages of gene summaries'.format(page_count))

            for index, (lsv_dict, genes) in enumerate(sg.get_paginated_genes_with_lsvs(args)):
                page_name = self.get_page_name(index)
                table_marks = tuple(table_marks_set(len(gene_set)) for gene_set in lsv_dict)
                next_page = self.get_next_page(index, page_count)
                group_names = metadata['group_names']

                self.add_to_voila_links(lsv_dict, page_name)
                log.debug('Write page {0}'.format(page_name))

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
                            lsvs=lsv_dict,
                            metadata=metadata,
                            database_name=database_name,
                            genome=genome
                    ):
                        html.write(el)

                prev_page = page_name

                log.debug('End summaries render')
