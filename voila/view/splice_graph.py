import json
import os
import uuid

import numpy

from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.exceptions import VoilaException
from voila.processes import VoilaPool
from voila.utils.voila_log import voila_log
from voila.view.html import Html, NumpyEncoder


class RenderSpliceGraphs(Html):

    def create_summary(self, paged):
        pass

    def __init__(self, args):
        super().__init__(args, None)
        self.db_id = uuid.uuid4().hex

        # self.view_metadata = m.view_metadata

        self.copy_static()
        self.render_sg_dbs()
        self.render_html(None, 'sg_summary.html')

    def sg_dbs(self, gene_ids):
        log = voila_log()
        for gene_id in gene_ids:
            with open(os.path.join(self.args.output, '{}.js'.format(gene_id)), 'w') as f:
                with ViewSpliceGraph(self.args) as sg:
                    metadata = sg.metadata
                    exp_name = metadata['experiment_names']

                    f.write('new PouchDB(\'voila_gene_{}\').bulkDocs(['.format(self.db_id))

                    log.debug('Write DB Gene ID: {}'.format(gene_id))

                    gene = sg.gene(gene_id)
                    gene_exp = sg.gene_experiment(gene, exp_name)
                    text = json.dumps(gene_exp)

                    f.write(text)
                    f.write(',')

                    del gene
                    del gene_exp
                    del text

                    f.write(']);')
                    f.write('\n')

    def render_sg_dbs(self):
        log = voila_log()

        log.debug('Create metadata file')
        with open(os.path.join(self.args.output, 'metadata.js'), 'w') as f:
            with ViewSpliceGraph(self.args) as sg:
                metadata = json.dumps(sg.metadata, cls=NumpyEncoder)
                f.write('new PouchDB(\'voila_gene_{}\').bulkDocs(['.format(self.db_id))
                f.write(metadata)
                f.write(',')
                f.write(']);')

                f.write('new PouchDB(\'voila_lsv_{}\').bulkDocs(['.format(self.db_id))
                f.write(metadata)
                f.write(',')
                f.write(']);')

                f.write('\n')

                f.write('const lsvs_arr = [')

                for gene_id in sg.gene_ids():
                    f.write(json.dumps({
                        'gene_id': gene_id
                    }))
                    f.write(',')
                f.write('];')

        with VoilaPool() as pool:
            # with self.matrix(self.args) as h:
            #     gene_ids = list(h.view_gene_ids())
            #     chunked_gene_ids = Html.chunkify(gene_ids, pool.processes)
            with ViewSpliceGraph(self.args) as sg:
                chunked_gene_ids = Html.chunkify(sg.gene_ids(), pool.processes)

            for p in [pool.apply_async(self.sg_dbs, (gene_ids,)) for gene_ids in chunked_gene_ids]:
                p.get()

    def fill_queue_gene_ids(self, queue, event):
        with ViewSpliceGraph(self.args) as sg:
            for gene in sg.genes():
                queue.put(gene.id)
        event.set()

    def db_genes(self, gene_ids, i):
        log = voila_log()
        try:
            args = self.args

            with ViewSpliceGraph(args) as sg, open(os.path.join(args.output, 'db_gene_{}.js'.format(i)),
                                                   'w') as db_gene:

                metadata = sg.metadata
                exp_name = metadata['experiment_names']

                db_gene.write('new PouchDB(\'voila_gene\').bulkDocs([')
                db_gene.write(json.dumps(metadata))
                db_gene.write(',')

                for gene_id in gene_ids:
                    log.debug('Write DB Gene ID: {}'.format(gene_id))

                    gene = sg.gene(gene_id)
                    gene_exp = sg.gene_experiment(gene, exp_name)
                    text = json.dumps(gene_exp)

                    db_gene.write(text)
                    db_gene.write(',')

                    del gene
                    del gene_exp
                    del text

                db_gene.write(']);')

        except Exception as e:
            log.exception(e)
            exit()

    # def render_dbs(self):
    #
    #     with VoilaPool() as pool:
    #         with ViewSpliceGraph(self.args) as sg:
    #             chunked_gene_ids = Html.chunkify([g.id for g in sg.genes()], pool.processes)
    #
    #         ps = [pool.apply_async(self.db_genes, (gene_ids, i)) for i, gene_ids in enumerate(chunked_gene_ids)]
    #         for p in ps:
    #             p.get()

    # def render_html(self):
    #     exec_dir = os.path.dirname(os.path.abspath(voila.__file__))
    #     template_dir = os.path.join(exec_dir, 'html')
    #     env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir), undefined=jinja2.StrictUndefined)
    #
    #     summary_template = env.get_template('sg_summary.html')
    #
    #     with open(os.path.join(self.args.output, 'summary.html'), 'w') as summary:
    #         summary.write(summary_template.render())

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
