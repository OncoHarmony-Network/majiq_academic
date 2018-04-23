import os

from voila.api import SpliceGraphs
from voila.utils.exceptions import VoilaException
from voila.utils.run_voila_utils import get_output_html, copy_static
from voila.utils.voila_log import voila_log
from voila.view.html import Html


class RenderSpliceGraphs(Html):
    def __init__(self, args):
        super(RenderSpliceGraphs, self).__init__(args)
        self.render_summaries()

    def render_summaries(self):
        voila_log().info('Rendering Splice Graph HTML output')
        summary_template = self.env.get_template('splice_graphs_summary_template.html')
        args = self.args
        output_html = get_output_html(args, args.splice_graph)
        summaries_subfolder = self.get_summaries_subfolder()
        log = voila_log()

        with SpliceGraphs(args.splice_graph, 'r') as sg:
            experiments = sg.get_experiments()
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
                            genes=genes,
                            experiments=experiments,
                            prev_page=prev_page,
                            next_page=next_page
                        )
                    )

                prev_page = page_name

        copy_static(args, index=False)
