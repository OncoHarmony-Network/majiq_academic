import os
from tempfile import NamedTemporaryFile

from voila import io_voila, constants
from voila.api import Voila, SpliceGraph
from voila.utils.exceptions import NoLsvsFound
from voila.utils.run_voila_utils import table_marks_set, copy_static, get_template_dir
from voila.utils.voila_log import voila_log
from voila.view.html import Html
from voila.voila_args import VoilaArgs


class Psi(Html, VoilaArgs):
    def __init__(self, args):
        """
        Render psi output.
        :param args: command line arguments
        :return: None
        """
        super(Psi, self).__init__(args)

        if not args.no_html:
            with Voila(args.voila_file, 'r') as v:
                self.metainfo = v.get_metainfo()
            self.render_summaries()
            self.render_index()
            copy_static(args)

        if not args.no_tsv:
            io_voila.tab_output(args, self.voila_links)

        if args.gtf:
            io_voila.generic_feature_format_txt_files(args)

        if args.gff:
            io_voila.generic_feature_format_txt_files(args, out_gff3=True)

    def render_index(self):
        log = voila_log()
        log.info('Render PSI HTML index')
        log.debug('Start index render')

        args = self.args
        env = self.env
        metainfo = self.metainfo
        index_row_template = env.get_template('psi_index_row.html')

        with NamedTemporaryFile(dir=get_template_dir()) as tmp_index_file:

            with Voila(args.voila_file, 'r') as v:

                lsv_count = v.get_lsv_count(args)
                too_many_lsvs = lsv_count > constants.MAX_LSVS_PSI_INDEX

                for index, lsv in enumerate(v.get_voila_lsvs(args)):
                    log.debug('Writing LSV {0} to index'.format(lsv.lsv_id))

                    for el in index_row_template.generate(
                            lsv=lsv,
                            index=index,
                            link=self.voila_links[lsv.name],
                            lexps=metainfo,
                            too_many_lsvs=too_many_lsvs
                    ):
                        tmp_index_file.write(bytearray(el, encoding='utf-8'))

                if not tmp_index_file.tell():
                    raise NoLsvsFound()

                with open(os.path.join(args.output, 'index.html'), 'w') as html:
                    index_template = env.get_template('index_single_summary_template.html')
                    for el in index_template.generate(
                            lexps=metainfo,
                            tmp_index_file_name=os.path.basename(tmp_index_file.name),
                            lsvs_count=lsv_count,
                            table_marks=table_marks_set(lsv_count),
                            prev_page=None,
                            next_page=None
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
        metainfo = self.metainfo
        group_name = metainfo['group_names'][0]

        with SpliceGraph(args.splice_graph, 'r') as sg:
            gene_experiments_list = sg.get_experiments()
            prev_page = None
            page_count = sg.get_page_count(args)

            for index, (lsv_dict, genes) in enumerate(sg.get_paginated_genes_with_lsvs(args)):
                page_name = self.get_page_name(index)
                table_marks = tuple(table_marks_set(len(gene_set)) for gene_set in lsv_dict)
                next_page = self.get_next_page(index, page_count)
                experiments = self.gene_experiments(metainfo['experiment_names'][0], genes, gene_experiments_list)

                self.add_to_voila_links(lsv_dict, page_name)

                log.debug('Writing {0}'.format(page_name))
                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:

                    for el in summary_template.generate(
                            table_marks=table_marks,
                            genes_dict=lsv_dict,
                            prev_page=prev_page,
                            next_page=next_page,
                            namePage=page_name,
                            lexps=metainfo,
                            genes_exps_list=experiments,
                            lsv_text_version=constants.LSV_TEXT_VERSION,
                            gtf=args.gtf,
                            group_name=group_name
                    ):
                        html.write(el)

                prev_page = page_name

        log.debug('End summaries render')
