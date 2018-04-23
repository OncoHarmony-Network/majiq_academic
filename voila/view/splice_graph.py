import os

import numpy

from voila.api.view_splice_graph import ViewSpliceGraph
from voila.exceptions import VoilaException
from voila.utils.voila_log import voila_log
from voila.view.html import Html


class RenderSpliceGraphs(Html):
    def render_dbs(self):
        pass

    def create_summary(self, paged):
        pass

    def __init__(self, args):
        super(RenderSpliceGraphs, self).__init__(args)
        self.copy_static(False)
        self.render_summaries()

    def render_summaries(self):
        voila_log().info('Rendering Splice Graph HTML output')
        summary_template = self.get_env().get_template('splice_graphs_summary_template.html')
        args = self.args
        output_html = self.get_output_html(args, args.splice_graph)
        summaries_subfolder = self.get_summaries_subfolder(args)
        log = voila_log()
        database_name = self.database_name()

        with ViewSpliceGraph(args) as sg:
            metadata = {'experiment_names': numpy.array([list(sg.experiment_names)]), 'group_names': [None]}
            prev_page = None
            page_count = sg.view_page_count()
            genome = sg.genome

            log.debug('Page count is {0}'.format(page_count))

            if not page_count:
                raise VoilaException('No Splice Graphs found')

            for index, genes in enumerate(sg.view_paginated_genes()):
                page_name = '{0}_{1}'.format(index, output_html)
                next_page = self.get_next_page(args, index, page_count)

                log.debug('Writing page {0}'.format(page_name))
                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:
                    html.write(
                        summary_template.render(
                            page_name=self.get_page_name(args, index),
                            genes=genes,
                            metadata=metadata,
                            prev_page=prev_page,
                            next_page=next_page,
                            database_name=database_name,
                            genome=genome
                        )
                    )

                prev_page = page_name
