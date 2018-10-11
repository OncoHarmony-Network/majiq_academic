import csv
import json
import os
from multiprocessing import Lock

import numpy as np

from voila import constants
from voila.api.view_matrix import ViewDeltaPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.exceptions import NotDeltaPsiVoilaFile
from voila.processes import VoilaPool
from voila.utils.voila_log import voila_log
from voila.view.html import Html, NumpyEncoder
from voila.view.tsv import Tsv
from voila.vlsv import matrix_area

lock = Lock()


class DeltaPsi(Html, Tsv):
    def __init__(self, args):
        super().__init__(args, ViewDeltaPsi)

        with ViewDeltaPsi(args) as m:
            if m.analysis_type != constants.ANALYSIS_DELTAPSI:
                raise NotDeltaPsiVoilaFile(args)
            self.view_metadata = m.view_metadata

        if not args.disable_html:
            self.copy_static()
            if not args.disable_db:
                self.render_dbs()
            self.render_html('dpsi_index.html', 'dpsi_summary.html')

        if not args.disable_tsv:
            self.delta_psi_tab_output()

    def dbs(self, gene_ids):
        log = voila_log()
        for gene_id in gene_ids:
            with open(os.path.join(self.args.output, '{}.js'.format(gene_id)), 'w') as f:
                with ViewDeltaPsi(self.args) as h, ViewSpliceGraph(self.args) as sg:
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

                            lsv = h.lsv(lsv_id).get_all()
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
            with ViewDeltaPsi(self.args) as h:
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
                    lsv = h.lsv(lsv_id)
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
            with ViewDeltaPsi(self.args) as h:
                gene_ids = list(h.view_gene_ids())
                chunked_gene_ids = Html.chunkify(gene_ids, pool.processes)

            for p in [pool.apply_async(self.dbs, (gene_ids,)) for gene_ids in chunked_gene_ids]:
                p.get()

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

                for gene_id in gene_ids:
                    gene = sg.gene(gene_id)
                    for lsv_id in m.view_gene_lsvs(gene.id):
                        log.debug('Write TSV row for {0}'.format(lsv_id))

                        lsv = m.lsv(lsv_id)
                        lsv_junctions = lsv.junctions
                        annot_juncs = sg.annotated_junctions(gene, lsv)
                        lsv_exons = sg.lsv_exons(gene, lsv)
                        excl_incl = list(lsv.excl_incl)
                        group_means = dict(lsv.group_means)

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
                            'De Novo Junctions': self.semicolon_join(annot_juncs),
                            'Junctions coords': self.semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'Exons coords': self.semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in self.filter_exons(lsv_exons)
                            ),
                            # TODO: this needs to be re-implemented
                            # 'IR coords': self.semicolon_join(
                            #     '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                            # ),
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
                                '%.3f' % i for i in group_means[group1]
                            ),
                            '%s E(PSI)' % group2: self.semicolon_join(
                                '%.3f' % i for i in group_means[group2]
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
        with ViewDeltaPsi():
            metadata = m.view_metadata

        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(dPSI) per LSV junction',
                      'P(|dPSI|>=%.2f) per LSV junction' % args.threshold,
                      'P(|dPSI|<=%.2f) per LSV junction' % args.non_changing_threshold,
                      '%s E(PSI)' % metadata['group_names'][0], '%s E(PSI)' % metadata['group_names'][1], 'LSV Type',
                      'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords']

        if voila_links:
            fieldnames.append('Voila link')

        self.write_tsv(fieldnames, ViewDeltaPsi)
