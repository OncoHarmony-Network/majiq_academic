import csv
import json
import os
from multiprocessing import Lock
from queue import Empty

import jinja2

import voila
from voila import constants
from voila.api.view_matrix import ViewHeterogens
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.exceptions import NotHeterogenVoilaFile
from voila.processes import VoilaPool, VoilaQueue
from voila.utils.voila_log import voila_log
from voila.view.html import Html, NumpyEncoder
from voila.view.tsv import Tsv
from voila.vlsv import get_expected_psi

tsv_lock = Lock()
gene_lock = Lock()
lsv_lock = Lock()


class Heterogen(Html, Tsv):
    def create_summary(self, paged):
        pass

    def __init__(self, args):
        super().__init__(args)
        with ViewHeterogens(args) as m:
            if m.analysis_type != constants.ANALYSIS_HETEROGEN:
                raise NotHeterogenVoilaFile(args)
            self.view_metadata = m.view_metadata

        if not args.disable_html:
            self.copy_static()
            self.render_html()
            if not args.disable_db:
                self.render_dbs()

        if not args.disable_tsv:
            self.render_tsv()

    def render_summaries(self):
        self.create_summaries(ViewHeterogens)

    def db_lsvs(self, q, e, i):
        log = voila_log()

        try:
            args = self.args
            with ViewHeterogens(args) as h:
                with open(os.path.join(args.output, 'db_lsv_{}.js'.format(i)), 'w') as db_lsv:

                    db_lsv.write('new PouchDB(\'voila_lsv\').bulkDocs([')

                    if i == 0:
                        db_lsv.write(json.dumps(h.view_metadata))
                        db_lsv.write(',')

                    while not (e.is_set() and q.empty()):
                        try:
                            lsv_id = q.get_nowait()
                        except Empty:
                            lsv_id = None

                        if lsv_id:
                            log.debug('Write DB LSV ID: {}'.format(lsv_id))
                            text = json.dumps(dict(h.heterogen(lsv_id).get_all()), cls=NumpyEncoder)
                            lsv_lock.acquire()
                            db_lsv.write(text)
                            db_lsv.write(',')
                            lsv_lock.release()
                            q.task_done()

                    db_lsv.write(']);')

        except Exception as e:
            log.exception(e)
            exit()

    def db_genes(self, gene_ids, i):
        log = voila_log()
        try:
            args = self.args

            with ViewHeterogens(args) as h:
                metadata = h.view_metadata

            with ViewSpliceGraph(args) as sg, open(os.path.join(args.output, 'db_gene_{}.js'.format(i)),
                                                   'w') as db_gene:

                exp_name = metadata['experiment_names']

                db_gene.write('new PouchDB(\'voila_gene\').bulkDocs([')

                if i == 0:
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

    def fill_queue_gene_ids(self, queue, event):
        with ViewHeterogens(self.args) as h:
            for gene_id in h.gene_ids:
                if any(h.view_gene_lsvs(gene_id)):
                    queue.put(gene_id)
        event.set()

    def fill_queue_lsv_ids(self, queue, event):
        event.clear()
        with ViewHeterogens(self.args) as h:
            for gene_id in h.gene_ids:
                for lsv_id in h.view_gene_lsvs(gene_id):
                    queue.put(lsv_id)
        event.set()

    def render_html(self):
        exec_dir = os.path.dirname(os.path.abspath(voila.__file__))
        template_dir = os.path.join(exec_dir, 'html')
        env = jinja2.Environment(extensions=["jinja2.ext.do"], loader=jinja2.FileSystemLoader(template_dir),
                                 undefined=jinja2.StrictUndefined)

        index_template = env.get_template('het_index.html')
        summary_template = env.get_template('het_summary.html')

        with open(os.path.join(self.args.output, 'index.html'), 'w') as index:
            index.write(index_template.render())

        with open(os.path.join(self.args.output, 'summary.html'), 'w') as summary:
            summary.write(summary_template.render())

    def render_dbs(self):
        args = self.args

        with VoilaPool() as pool:
            with ViewSpliceGraph(self.args) as sg:
                chunked_gene_ids = Html.chunkify([g.id for g in sg.genes()], pool.processes)

            ps = [pool.apply_async(self.db_genes, (gene_ids, i)) for i, gene_ids in enumerate(chunked_gene_ids)]
            for p in ps:
                p.get()

            with VoilaQueue(self.fill_queue_lsv_ids) as (q, e):
                for p in [pool.apply_async(self.db_lsvs, (q, e, i)) for i in range(pool.processes)]:
                    p.get()

        with open(os.path.join(args.output, 'db_lsv.js'), 'a') as db_lsv:
            db_lsv.write('])};')

    def tsv_row(self, q, e, tsv_file, fieldnames):
        args = self.args
        log = voila_log()
        metadata = self.view_metadata
        group_names = metadata['group_names']
        with ViewHeterogens(args) as m, ViewSpliceGraph(args) as sg:
            with open(tsv_file, 'a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                while not (e.is_set() and q.empty()):
                    try:
                        lsv_id = q.get_nowait()
                    except Empty:
                        lsv_id = None

                    if lsv_id:
                        lsv = m.heterogen(lsv_id)
                        gene = sg.gene(lsv.gene_id)
                        log.debug('Write TSV row for {0}'.format(lsv_id))
                        lsv = m.heterogen(lsv_id)
                        lsv_junctions = list(gene.lsv_junctions(lsv))
                        lsv_exons = list(gene.lsv_exons(lsv))
                        mean_psi = list(lsv.mean_psi)

                        row = {
                            'Gene Name': gene.name,
                            'Gene ID': gene.id,
                            'LSV ID': lsv_id,
                            'LSV Type': lsv.lsv_type,
                            'A5SS': lsv.a5ss,
                            'A3SS': lsv.a3ss,
                            'ES': lsv.exon_skipping,
                            'Num. Junctions': len(lsv_junctions),
                            'Num. Exons': lsv.exon_count,
                            'chr': gene.chromosome,
                            'strand': gene.strand,
                            'De Novo Junctions': self.semicolon_join(
                                int(not junc.annotated) for junc in lsv_junctions
                            ),
                            'Junctions coords': self.semicolon_join(
                                '{0}-{1}'.format(junc.start, junc.end) for junc in lsv_junctions
                            ),
                            'Exons coords': self.semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in self.filter_exons(lsv_exons)
                            ),
                            # todo: reimplement this column
                            # 'IR coords': self.semicolon_join(
                            #     '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                            # ),
                        }

                        for idx, group in enumerate(group_names):
                            if mean_psi[idx] is not None:
                                row['%s E(PSI)' % group] = self.semicolon_join(
                                    get_expected_psi(x) for x in mean_psi[idx])

                        for key, value in lsv.junction_stats:
                            row[key] = self.semicolon_join(value)

                        # if voila_links:
                        #     summary_path = voila_links[gene_id]
                        #     if not os.path.isabs(summary_path):
                        #         summary_path = join(os.getcwd(), args.output, summary_path)
                        #     row['Voila link'] = "file://{0}".format(summary_path)

                        tsv_lock.acquire()
                        writer.writerow(row)
                        tsv_lock.release()
                        q.task_done()

    def render_tsv(self):
        log = voila_log()

        with ViewHeterogens(self.args) as m:
            metadata = m.view_metadata
            fieldnames = ['Gene Name', 'Gene ID', 'LSV ID', 'LSV Type', 'strand', 'chr'] + \
                         ['%s E(PSI)' % group for group in metadata['group_names']] + \
                         list(m.junction_stats_column_names) + \
                         ['A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions',
                          'Junctions coords', 'Exons coords', 'IR coords']

        log.info("Creating Tab-delimited output file")

        args = self.args
        output_html = Html.get_output_html(args, args.voila_file[0])
        tsv_file = os.path.join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

        with open(tsv_file, 'w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

        with VoilaPool() as vp, VoilaQueue(self.fill_queue_lsv_ids) as (q, e):
            for x in [vp.apply_async(self.tsv_row, (q, e, tsv_file, fieldnames)) for _ in range(vp.processes)]:
                x.get()

        log.info("Delimited output file successfully created in: %s" % tsv_file)
