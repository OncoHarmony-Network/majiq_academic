import csv
import errno
import os
from multiprocessing import Lock

import numpy as np

from voila import constants
from voila.api.matrix_hdf5 import lsv_id_to_gene_id
from voila.api.view_matrix import ViewDeltaPsi
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.exceptions import NotDeltaPsiVoilaFile
from voila.utils.voila_log import voila_log
from voila.utils.voila_pool import VoilaPool
from voila.view.html import Html
from voila.view.tsv import Tsv
from voila.vlsv import matrix_area

lock = Lock()


class DeltaPsi(Html, Tsv):
    def __init__(self, args):
        super(DeltaPsi, self).__init__(args)
        with ViewDeltaPsi(args) as m:
            if m.analysis_type != constants.ANALYSIS_DELTAPSI:
                raise NotDeltaPsiVoilaFile(args.voila_file)
            self.view_metadata = m.view_metadata

        if not args.disable_html:
            self.copy_static(args)
            self.create_db_files()
            self.render_summaries()
            self.render_index()

        if not args.disable_tsv:
            self.delta_psi_tab_output()

    def render_index(self):
        log = voila_log()
        log.info('Render Delta PSI HTML index')
        log.debug('Start index render')
        args = self.args
        metadata = self.view_metadata

        with ViewSpliceGraph(args) as sg, ViewDeltaPsi(args) as m:
            lsv_count = m.view_lsv_count()
            too_many_lsvs = lsv_count > constants.MAX_LSVS_DELTAPSI_INDEX

            with open(os.path.join(args.output, 'index.html'), 'w') as html:
                index_template = self.get_env().get_template('index_delta_summary_template.html')
                html.write(
                    index_template.render(
                        lexps=metadata,
                        table_marks=self.table_marks_set(lsv_count),
                        lsvs_count=lsv_count,
                        prev_page=None,
                        next_page=None,
                        database_name='deltapsi_' + self.database_name(),
                        lsvs=list(m.delta_psi(lsv_id) for lsv_id in m.view_lsv_ids()),
                        genes=sg.gene,
                        gene_ids=set(lsv_id_to_gene_id(lsv_id) for lsv_id in m.view_lsv_ids()),
                        links=self.voila_links,
                        too_many_lsvs=too_many_lsvs,
                        metadata=metadata,
                        genome=sg.genome,
                        lsv_text_version=constants.LSV_TEXT_VERSION,
                        threshold=args.threshold
                    )
                )

                log.debug('End index render')

    @classmethod
    def create_gene_db(cls, gene_ids, args, experiment_names):
        template = cls.get_env().get_template('gene_db_template.html')
        log = voila_log()
        with ViewSpliceGraph(args) as sg, ViewDeltaPsi(args) as m:
            for gene in sg.genes(gene_ids):
                log.debug('creating {}'.format(gene.id))

                with open(os.path.join(args.output, 'db', '{}.js'.format(gene.id)), 'w') as html:
                    html.write(
                        template.render(
                            gene=gene.get_experiment(experiment_names),
                            lsvs=tuple(m.delta_psi(lsv_id) for lsv_id in m.view_gene_lsvs(gene.id))
                        )
                    )

    @classmethod
    def create_summary(cls, metadata, args, database_name, paged):
        summary_template = cls.get_env().get_template("deltapsi_summary_template.html")
        summaries_subfolder = cls.get_summaries_subfolder(args)
        group_names = metadata['group_names']
        links = {}

        with ViewSpliceGraph(args) as sg, ViewDeltaPsi(args) as m:
            genome = sg.genome
            page_count = m.page_count()

            for index, gene_ids in paged:
                page_name = cls.get_page_name(args, index)
                next_page = cls.get_next_page(args, index, page_count)
                prev_page = cls.get_prev_page(args, index)

                lsv_dict = {gene_id: tuple(lsv_id for lsv_id in m.view_gene_lsvs(gene_id)) for
                            gene_id in gene_ids}
                table_marks = tuple(cls.table_marks_set(len(gene_set)) for gene_set in lsv_dict)

                with open(os.path.join(summaries_subfolder, page_name), 'w') as html:
                    html.write(
                        summary_template.render(
                            page_name=page_name,
                            threshold=args.threshold,
                            lsv_text_version=constants.LSV_TEXT_VERSION,
                            table_marks=table_marks,
                            prev_page=prev_page,
                            next_page=next_page,
                            gtf=args.gtf,
                            group_names=group_names,
                            genes=list(sg.genes(gene_ids)),
                            lsv_ids=lsv_dict,
                            delta_psi_lsv=m.delta_psi,
                            metadata=metadata,
                            database_name=database_name,
                            genome=genome,
                            splice_graph=sg
                        )
                    )

                links.update(dict(cls.voila_links(lsv_dict, page_name)))

            return links

    def render_summaries(self):
        log = voila_log()
        log.info('Render Delta PSI HTML summaries')
        log.debug('Start summaries render')

        args = self.args
        metadata = self.view_metadata
        database_name = self.database_name()

        with ViewDeltaPsi(args) as m:
            paged_genes = tuple(m.paginated_genes())

        multiple_results = []
        with VoilaPool() as vp:
            for paged in self.chunkify(tuple(enumerate(paged_genes)), vp.processes):
                multiple_results.append(
                    vp.pool.apply_async(self.create_summary, (metadata, args, database_name, paged)))

            for res in multiple_results:
                self.voila_links.update(res.get())

    def create_db_files(self):
        args = self.args
        metadata = self.view_metadata
        log = voila_log()
        log.info('Create DB files')

        with ViewDeltaPsi(args) as m:
            gene_ids = tuple(m.view_gene_ids())

        try:
            os.makedirs(os.path.join(args.output, 'db'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        def chunkify(lst, n):
            for i in range(n):
                yield lst[i::n]

        names = metadata['experiment_names']
        multiple_results = []

        with VoilaPool() as vp:
            for genes in chunkify(gene_ids, vp.processes):
                multiple_results.append(vp.apply_async(self.create_gene_db, (genes, args, names)))

            for res in multiple_results:
                res.get()

        log.debug('finished writing db files.')

    def tsv_row(self, gene_ids, tsv_file, fieldnames):
        voila_links = self.voila_links
        args = self.args
        log = voila_log()

        with ViewDeltaPsi(args) as m, ViewSpliceGraph(args) as sg:
            metadata = m.view_metadata
            group1 = metadata['group_names'][0]
            group2 = metadata['group_names'][1]

            with open(tsv_file, 'a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene in sg.genes(gene_ids):
                    for lsv_id in m.view_gene_lsvs(gene.id):
                        log.debug('Write TSV row for {0}'.format(lsv_id))
                        lsv = m.delta_psi(lsv_id)

                        lsv_junctions = list(gene.lsv_junctions(lsv))
                        lsv_exons = list(gene.lsv_exons(lsv))
                        group_means = list(lsv.group_means)
                        excl_incl = list(lsv.excl_incl)

                        row = {
                            '#Gene Name': gene.name,
                            'Gene ID': gene.id,
                            'LSV ID': lsv_id,
                            'LSV Type': lsv.lsv_type,
                            'A5SS': lsv.prime5,
                            'A3SS': lsv.prime3,
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
                            'IR coords': self.semicolon_join(
                                '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                            ),
                            'E(dPSI) per LSV junction': self.semicolon_join(
                                excl_incl[i][1] - excl_incl[i][0] for i in
                                range(np.size(lsv.bins, 0))
                            ),
                            'P(|dPSI|>=%.2f) per LSV junction' % args.threshold: self.semicolon_join(
                                matrix_area(b, args.threshold) for b in lsv.bins
                            ),
                            'P(|dPSI|<=%.2f) per LSV junction' % args.non_changing_threshold: self.semicolon_join(
                                lsv.high_probability_non_changing()
                            ),
                            '%s E(PSI)' % group1: self.semicolon_join(
                                '%.3f' % i for i in group_means[0]
                            ),
                            '%s E(PSI)' % group2: self.semicolon_join(
                                '%.3f' % i for i in group_means[1]
                            )
                        }

                        if voila_links:
                            summary_path = voila_links[gene.id]
                            if not os.path.isabs(summary_path):
                                summary_path = os.path.join(os.getcwd(), args.output, summary_path)
                            row['Voila link'] = "file://{0}".format(summary_path)

                        lock.acquire()
                        writer.writerow(row)
                        lock.release()

    def delta_psi_tab_output(self):
        args = self.args
        voila_links = self.voila_links
        metadata = self.view_metadata

        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(dPSI) per LSV junction',
                      'P(|dPSI|>=%.2f) per LSV junction' % args.threshold,
                      'P(|dPSI|<=%.2f) per LSV junction' % args.non_changing_threshold,
                      '%s E(PSI)' % metadata['group_names'][0], '%s E(PSI)' % metadata['group_names'][1], 'LSV Type',
                      'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords']

        if voila_links:
            fieldnames.append('Voila link')

        self.write_tsv(fieldnames)
