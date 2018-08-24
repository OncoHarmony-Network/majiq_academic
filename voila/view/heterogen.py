import csv
import glob
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
meta_lock = Lock()


class Heterogen(Html, Tsv):
    def create_summary(self, paged):
        pass

    def __init__(self, args):
        super().__init__(args, ViewHeterogens)
        with ViewHeterogens(args) as m:
            print(m.analysis_type)
            if m.analysis_type != constants.ANALYSIS_HETEROGEN:
                raise NotHeterogenVoilaFile(args)
            self.view_metadata = m.view_metadata

        if not args.disable_html:
            self.copy_static()
            if not args.disable_db:
                self.render_dbs()
            self.render_html()

        if not args.disable_tsv:
            self.render_tsv()

    def render_summaries(self):
        self.create_summaries(ViewHeterogens)

    def dbs(self, gene_ids):
        log = voila_log()
        for gene_id in gene_ids:
            with open(os.path.join(self.args.output, '{}.js'.format(gene_id)), 'w') as f:
                with ViewHeterogens(self.args) as h, ViewSpliceGraph(self.args) as sg:
                    metadata = h.view_metadata
                    exp_name = metadata['experiment_names']
                    lsv_ids = h.view_gene_lsvs(gene_id)

                    f.write('new PouchDB(\'voila_gene\').bulkDocs([')

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

                    if lsv_ids:
                        f.write('new PouchDB(\'voila_lsv\').bulkDocs([')

                        for lsv_id in lsv_ids:
                            log.debug('Write DB LSV ID: {}'.format(lsv_id))

                            lsv = h.heterogen(lsv_id).get_all()
                            lsv_dict = dict(lsv)
                            text = json.dumps(lsv_dict, cls=NumpyEncoder)

                            f.write(text)
                            f.write(',')

                            del lsv
                            del lsv_dict
                            del text

                        f.write(']);')

    def db_lsvs(self, lsv_ids, i):
        log = voila_log()

        try:
            lsv_ids_list = self.chunks(lsv_ids, 50)
            args = self.args

            for page_idx, lsv_ids in enumerate(lsv_ids_list):

                with ViewHeterogens(args) as h:
                    with open(os.path.join(args.output, 'db_lsv_{}_{}.js'.format(i, page_idx)), 'w') as db_lsv:

                        db_lsv.write('new PouchDB(\'voila_lsv\').bulkDocs([')

                        if i == 0 and page_idx == 0:
                            db_lsv.write(json.dumps(h.view_metadata))
                            db_lsv.write(',')

                        for lsv_id in lsv_ids:
                            log.debug('Write DB LSV ID: {}'.format(lsv_id))

                            lsv = h.heterogen(lsv_id).get_all()
                            lsv_dict = dict(lsv)
                            text = json.dumps(lsv_dict, cls=NumpyEncoder)

                            db_lsv.write(text)
                            db_lsv.write(',')

                        db_lsv.write(']);')

        except Exception as e:
            log.exception(e)
            exit()

    @staticmethod
    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def db_genes(self, gene_ids, i):
        log = voila_log()

        try:
            gene_ids_list = self.chunks(gene_ids, 50)
            args = self.args

            with ViewHeterogens(args) as h:
                metadata = h.view_metadata
                exp_name = metadata['experiment_names']

            for page_idx, gene_ids in enumerate(gene_ids_list):

                with ViewSpliceGraph(args) as sg:
                    with open(os.path.join(args.output, 'db_gene_{}_{}.js'.format(i, page_idx)), 'w') as db_gene:

                        db_gene.write('new PouchDB(\'voila_gene\').bulkDocs([')

                        if i == 0 and page_idx == 0:
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

        db_lsvs = [os.path.basename(f) for f in glob.glob(os.path.join(self.args.output, 'db_lsv*.js'))]
        db_genes = [os.path.basename(f) for f in glob.glob(os.path.join(self.args.output, 'db_gene*.js'))]
        dbs = [os.path.basename(f) for f in glob.glob(os.path.join(self.args.output, '*.js'))]

        with open(os.path.join(self.args.output, 'index.html'), 'w') as index:
            index.write(index_template.render({
                'db_lsvs': db_lsvs,
                'db_genes': db_genes,
                'dbs': dbs
            }))

        with open(os.path.join(self.args.output, 'summary.html'), 'w') as summary:
            summary.write(summary_template.render({
                'db_lsvs': db_lsvs,
                'db_genes': db_genes
            }))

    def render_dbs(self):
        log = voila_log()

        log.debug('Create metadata file')
        with open(os.path.join(self.args.output, 'metadata.js'), 'w') as f:
            with ViewHeterogens(self.args) as h:
                metadata = json.dumps(h.view_metadata)
                f.write('new PouchDB(\'voila_gene\').bulkDocs([')
                f.write(metadata)
                f.write(',')
                f.write(']);')

                f.write('new PouchDB(\'voila_lsv\').bulkDocs([')
                f.write(metadata)
                f.write(',')
                f.write(']);')

                f.write('\n')

                f.write('const lsvs_arr = [')

                for lsv_id in h.view_lsv_ids():
                    lsv = h.heterogen(lsv_id)
                    f.write(json.dumps({
                        '_id': lsv.lsv_id,
                        'target': lsv.target,
                        'binary': lsv.binary,
                        'exon_skipping': lsv.exon_skipping,
                        'A5SS': lsv.a5ss,
                        'A3SS': lsv.a3ss,
                        'gene_id': lsv.gene_id
                    }))
                    f.write(',')
                f.write('];')

        with VoilaPool() as pool:
            with ViewHeterogens(self.args) as h:
                gene_ids = list(h.view_gene_ids())
                chunked_gene_ids = Html.chunkify(gene_ids, pool.processes)

            for p in [pool.apply_async(self.dbs, (gene_ids,)) for gene_ids in chunked_gene_ids]:
                p.get()

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
