import os

import numpy

from voila.api.view_splice_graph import ViewSpliceGraph
from voila.utils.exceptions import VoilaException
from voila.utils.run_voila_utils import get_output_html, copy_static
from voila.utils.voila_log import voila_log
from voila.view.html import Html
from voila.voila_args import VoilaArgs


class RenderSpliceGraphs(Html, VoilaArgs):
    def __init__(self, args):
        super(RenderSpliceGraphs, self).__init__(args)
        self.render_summaries()
        copy_static(args, index=False)

    @classmethod
    def arg_parents(cls):
        # base, gene_search, parser_splicegraphs, output

        parser = cls.get_parser()

        parser.add_argument('splice_graph',
                            type=cls.check_splice_graph_file,
                            help=cls.SPLICE_GRAPH_HELP)

        parser.add_argument('--limit',
                            type=int,
                            default=20,
                            help='Limit the number of splice graphs shown.  Default is 20.')

        return (
            cls.base_args(), cls.gene_search_args(), cls.output_args(), parser
        )

    def render_summaries(self):
        voila_log().info('Rendering Splice Graph HTML output')
        summary_template = self.env.get_template('splice_graphs_summary_template.html')
        args = self.args
        output_html = get_output_html(args, args.splice_graph)
        summaries_subfolder = self.get_summaries_subfolder()
        log = voila_log()
        database_name = self.database_name()

        with ViewSpliceGraph(args.splice_graph, 'r') as sg:
            metadata = {'experiment_names': numpy.array([list(sg.experiment_names)]), 'group_names': [None]}
            prev_page = None
            page_count = sg.get_page_count(args)

            log.debug('Page count is {0}'.format(page_count))

            if not page_count:
                raise VoilaException('No Splice Graphs found')

            for index, genes in enumerate(sg.get_paginated_genes(args)):
                page_name = '{0}_{1}'.format(index, output_html)
                next_page = self.get_next_page(index, page_count)

                log.debug('Writing page {0}'.format(page_name))
                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:
                    html.write(
                        summary_template.render(
                            page_name=self.get_page_name(index),
                            genes=[sg.gene(gene_id) for gene_id in genes],
                            metadata=metadata,
                            prev_page=prev_page,
                            next_page=next_page,
                            database_name=database_name,
                        )
                    )

                prev_page = page_name
