import json
import multiprocessing
import os
import queue

import jinja2
import numpy

import voila
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.exceptions import VoilaException
from voila.processes import VoilaPool, VoilaQueue
from voila.utils.voila_log import voila_log
from voila.view.html import Html

gene_lock = multiprocessing.Lock()


class RenderSpliceGraphs(Html):

    def create_summary(self, paged):
        pass

    def __init__(self, args):
        super(RenderSpliceGraphs, self).__init__(args)
        self.copy_static()
        self.render_html()
        if not args.disable_db:
            self.render_dbs()

    def fill_queue_gene_ids(self, queue, event):
        with ViewSpliceGraph(self.args) as h:
            for gene_id in h.gene_ids:
                queue.put(gene_id)
        event.set()

    def db_genes(self, q, e):
        log = voila_log()
        try:
            args = self.args

            with ViewSpliceGraph(args) as sg:

                metadata = sg.metadata

                with open(os.path.join(args.output, 'db_gene.js'), 'a') as db_gene:
                    while not (e.is_set() and q.empty()):
                        try:
                            gene_id = q.get_nowait()
                        except queue.Empty:
                            gene_id = None

                        if gene_id:
                            log.debug('Write DB Gene ID: {}'.format(gene_id))
                            text = json.dumps(sg.gene(gene_id).get_experiment(metadata['experiment_names']))
                            gene_lock.acquire()
                            db_gene.write(text)
                            db_gene.write(',')
                            gene_lock.release()

                            q.task_done()
        except Exception as e:
            log.exception(e)
            exit()

    def render_dbs(self):
        args = self.args

        with ViewSpliceGraph(args) as sg:
            metadata = sg.metadata

        with open(os.path.join(args.output, 'db_gene.js'), 'w') as db_gene:
            db_gene.write('db_gene.bulkDocs([')
            db_gene.write(json.dumps(metadata))
            db_gene.write(',')

        with VoilaPool() as pool:
            with VoilaQueue(self.fill_queue_gene_ids) as (q, e):
                for p in [pool.apply_async(self.db_genes, (q, e)) for _ in range(pool.processes)]:
                    p.get()

        with open(os.path.join(args.output, 'db_gene.js'), 'a') as db_gene:
            db_gene.write('])')

    def render_html(self):
        exec_dir = os.path.dirname(os.path.abspath(voila.__file__))
        template_dir = os.path.join(exec_dir, 'html')
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir), undefined=jinja2.StrictUndefined)

        summary_template = env.get_template('sg_summary.html')

        with open(os.path.join(self.args.output, 'summary.html'), 'w') as summary:
            summary.write(summary_template.render())

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
