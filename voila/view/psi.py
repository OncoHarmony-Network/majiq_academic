import csv
import glob
import json
import os
import uuid
from multiprocessing import Lock

import jinja2

import voila
from voila import constants
from voila.api.matrix_hdf5 import lsv_id_to_gene_id
from voila.api.view_matrix import ViewPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.exceptions import NotPsiVoilaFile
from voila.processes import VoilaPool
from voila.utils.voila_log import voila_log
from voila.view.html import Html, NumpyEncoder
from voila.view.tsv import Tsv

lock = Lock()


class Psi(Html, Tsv):
    def __init__(self, args):
        super().__init__(args)
        self.db_id = uuid.uuid4().hex

        with ViewPsi(args) as m:
            if m.analysis_type != constants.ANALYSIS_PSI:
                raise NotPsiVoilaFile(args)
            self.view_metadata = m.view_metadata

        if not args.disable_html:
            self.copy_static()
            if not args.disable_db:
                self.render_dbs()
            self.render_html()

        if not args.disable_tsv:
            self.render_tsv()

    def render_index(self):
        log = voila_log()
        log.info('Render PSI HTML index')
        log.debug('Start index render')

        args = self.args
        metadata = self.view_metadata

        with ViewPsi(args) as m, ViewSpliceGraph(args) as sg:
            lsv_ids = tuple(m.view_lsv_ids())
            lsv_count = m.view_lsv_count()
            too_many_lsvs = lsv_count > constants.MAX_LSVS_PSI_INDEX

            with open(os.path.join(args.output, 'het_index.html'), 'w') as html:
                index_template = self.get_env().get_template('index_single_summary_template.html')
                html.write(
                    index_template.render(
                        lexps=metadata,
                        metadata=metadata,
                        lsvs_count=lsv_count,
                        table_marks=self.table_marks_set(lsv_count),
                        prev_page=None,
                        next_page=None,
                        links=self.voila_links,
                        genes=sg.gene,
                        gene_ids=set(lsv_id_to_gene_id(lsv_id) for lsv_id in m.view_lsv_ids()),
                        lsvs=(m.psi(lsv_id) for lsv_id in lsv_ids),
                        too_many_lsvs=too_many_lsvs,
                        database_name=self.database_name(),
                        genome=sg.genome,
                        lsv_text_version=constants.LSV_TEXT_VERSION
                    )
                )

    def create_summary(self, paged):
        args = self.args
        summary_template = self.get_env().get_template("psi_summary_template.html")
        summaries_subfolder = self.get_summaries_subfolder(args)
        database_name = self.database_name()
        metadata = self.view_metadata
        links = {}

        with ViewSpliceGraph(args) as sg, ViewPsi(args) as m:
            genome = sg.genome
            page_count = m.page_count()

            for index, genes in paged:
                page_name = self.get_page_name(args, index)
                next_page = self.get_next_page(args, index, page_count)
                prev_page = self.get_prev_page(args, index)
                lsv_dict = {gene_id: tuple(lsv_id for lsv_id in m.view_gene_lsvs(gene_id)) for
                            gene_id in genes}
                table_marks = tuple(self.table_marks_set(len(gene_set)) for gene_set in lsv_dict)

                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:
                    html.write(
                        summary_template.render(
                            genes=list(sg.gene(gene_id) for gene_id in genes),
                            table_marks=table_marks,
                            lsv_ids=lsv_dict,
                            psi_lsv=m.psi,
                            prev_page=prev_page,
                            next_page=next_page,
                            namePage=page_name,
                            metadata=metadata,
                            lsv_text_version=constants.LSV_TEXT_VERSION,
                            gtf=args.gtf,
                            database_name=database_name,
                            genome=genome
                        )
                    )

                links.update(self.get_voila_links(lsv_dict, page_name))

            return links

    def render_html(self):
        exec_dir = os.path.dirname(os.path.abspath(voila.__file__))
        template_dir = os.path.join(exec_dir, 'html')
        env = jinja2.Environment(extensions=["jinja2.ext.do"], loader=jinja2.FileSystemLoader(template_dir),
                                 undefined=jinja2.StrictUndefined)

        index_template = env.get_template('psi_index.html')
        summary_template = env.get_template('psi_summary.html')

        with open(os.path.join(self.args.output, 'index.html'), 'w') as index:
            index.write(index_template.render({
                'db_id': self.db_id
            }))

        with open(os.path.join(self.args.output, 'summary.html'), 'w') as summary:
            summary.write(summary_template.render({
                'db_id': self.db_id
            }))

    def dbs(self, gene_ids):
        log = voila_log()
        for gene_id in gene_ids:
            with open(os.path.join(self.args.output, '{}.js'.format(gene_id)), 'w') as f:
                with ViewPsi(self.args) as h, ViewSpliceGraph(self.args) as sg:
                    metadata = h.view_metadata
                    exp_name = metadata['experiment_names']
                    lsv_ids = h.view_gene_lsvs(gene_id)

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

                    if lsv_ids:
                        f.write('new PouchDB(\'voila_lsv_{}\').bulkDocs(['.format(self.db_id))

                        for lsv_id in lsv_ids:
                            log.debug('Write DB LSV ID: {}'.format(lsv_id))

                            lsv = h.psi(lsv_id).get_all()
                            lsv_dict = dict(lsv)
                            text = json.dumps(lsv_dict, cls=NumpyEncoder)

                            f.write(text)
                            f.write(',')

                            del lsv
                            del lsv_dict
                            del text

                        f.write(']);')

    def render_dbs(self):
        log = voila_log()

        log.debug('Create metadata file')
        with open(os.path.join(self.args.output, 'metadata.js'), 'w') as f:
            with ViewPsi(self.args) as h:
                metadata = json.dumps(h.view_metadata, cls=NumpyEncoder)
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

                for lsv_id in h.view_lsv_ids():
                    lsv = h.psi(lsv_id)
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
            with ViewPsi(self.args) as h:
                gene_ids = list(h.view_gene_ids())
                chunked_gene_ids = Html.chunkify(gene_ids, pool.processes)

            for p in [pool.apply_async(self.dbs, (gene_ids,)) for gene_ids in chunked_gene_ids]:
                p.get()

    def tsv_row(self, gene_ids, tsv_file, fieldnames):
        args = self.args
        voila_links = self.voila_links
        log = voila_log()
        with ViewPsi(args) as m, ViewSpliceGraph(args) as sg:

            with open(tsv_file, 'a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene in sg.genes(gene_ids):
                    for lsv_id in m.view_gene_lsvs(gene.id):
                        lsv = m.psi(lsv_id)
                        lsv_junctions = list(gene.lsv_junctions(lsv))
                        lsv_exons = list(gene.lsv_exons(lsv))

                        row = {
                            '#Gene Name': gene.name,
                            'Gene ID': gene.id,
                            'LSV ID': lsv_id,
                            'LSV Type': lsv.lsv_type,
                            'A5SS': lsv.a5ss,
                            'A3SS': lsv.a3ss,
                            'ES': lsv.exon_skipping,
                            'Num. Junctions': lsv.junction_count,
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
                            # 'IR coords': self.semicolon_join(
                            #     '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                            # ),
                            'E(PSI) per LSV junction': self.semicolon_join(lsv.means),
                            'Var(E(PSI)) per LSV junction': self.semicolon_join(lsv.variances)
                        }

                        if voila_links:
                            summary_path = voila_links[gene.id]
                            if not os.path.isabs(summary_path):
                                summary_path = os.path.join(os.getcwd(), args.output, summary_path)
                            row['Voila link'] = "file://{0}".format(summary_path)

                        log.debug('Write TSV row for {0}'.format(lsv_id))

                        lock.acquire()
                        writer.writerow(row)
                        lock.release()

    def psi_tab_output(self):
        voila_links = self.voila_links

        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction',
                      'LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords']
        if voila_links:
            fieldnames.append('Voila link')

        self.write_tsv(fieldnames, ViewPsi)
