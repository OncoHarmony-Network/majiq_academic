import csv
import multiprocessing
import os
from _queue import Empty
from abc import ABC

import numpy as np

from voila import constants
from voila.api.view_matrix import ViewPsi, ViewDeltaPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.config import Config
from voila.processes import VoilaPool, VoilaQueue
from voila.utils.voila_log import voila_log
from voila.view.html import Html
from voila.vlsv import matrix_area

lock = multiprocessing.Lock()


class Tsv(ABC):
    def __init__(self, args):
        self.args = args

    @staticmethod
    def semicolon_join(value_list):
        return ';'.join(str(x) for x in value_list)

    @staticmethod
    def filter_exons(exons):
        for exon in exons:
            if exon[0] == -1:
                yield 'nan', exon[1]
            elif exon[1] == -1:
                yield exon[0], 'nan'
            else:
                yield exon

    def write_tsv(self, fieldnames, view_matrix):
        log = voila_log()
        log.info("Creating Tab-delimited output file")

        args = self.args
        output_html = Html.get_output_html(args, args.voila_files[0])
        tsv_file = os.path.join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

        with view_matrix(args) as m:
            view_gene_ids = list(m.gene_ids)

        with open(tsv_file, 'w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

        multiple_results = []

        with VoilaPool() as vp, VoilaQueue() as vq:
            for gene_ids in self.chunkify(view_gene_ids, vp.processes):
                multiple_results.append(vp.apply_async(self.tsv_row, (gene_ids, tsv_file, fieldnames)))

        [r.get() for r in multiple_results]

        log.info("Delimited output file successfully created in: %s" % tsv_file)


def chunkify(lst, n):
    for i in range(n):
        yield lst[i::n]


def filter_exons(exons):
    for exon in exons:
        if exon[0] == -1:
            yield 'nan', exon[1]
        elif exon[1] == -1:
            yield exon[0], 'nan'
        else:
            yield exon


def semicolon_join(value_list):
    return ';'.join(str(x) for x in value_list)


class NewTsv:
    def __init__(self):
        self.config = Config()
        analysis_type = self.config.analysis_type
        if analysis_type == constants.ANALYSIS_PSI:
            PsiTsv()
        elif analysis_type == constants.ANALYSIS_DELTAPSI:
            DeltaPsiTsv()


class AnalysisTypeTsv:
    def __init__(self, view_matrix):
        self.config = Config()
        self.view_matrix = view_matrix

    @staticmethod
    def tsv_row(self, q, e, tsv_file, fieldnames):
        raise NotImplementedError()

    def fill_queue(self, queue, event):
        with self.view_matrix() as m:
            for gene_id in m.gene_ids:
                # # if any(h.view_gene_lsvs(gene_id)):
                # #     queue.put(gene_id)
                queue.put(gene_id)
            event.set()

    def gene_ids(self, q, e):
        while not (e.is_set() and q.empty()):
            try:
                yield q.get_nowait()
            except Empty:
                pass

        assert q.empty()

    def write_tsv(self, fieldnames):
        config = self.config
        log = voila_log()
        log.info("Creating Tab-delimited output file")
        nproc = self.config.nproc
        multiple_results = []
        output_html = Html.get_output_html(config.output, config.voila_file)
        tsv_file = os.path.join(config.output, output_html.rsplit('.html', 1)[0] + '.tsv')
        mgr = multiprocessing.Manager()
        queue = mgr.Queue(nproc * 2)
        event = mgr.Event()
        fill_queue_proc = multiprocessing.Process(target=self.fill_queue, args=(queue, event))

        log.debug('Manager PID {}'.format(mgr._process.ident))

        fill_queue_proc.start()

        with open(tsv_file, 'w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

        with multiprocessing.Pool(processes=nproc) as pool:
            for _ in range(nproc):
                multiple_results.append(pool.apply_async(self.tsv_row, (queue, event, tsv_file, fieldnames)))
            [r.get() for r in multiple_results]

        log.info("Delimited output file successfully created in: %s" % tsv_file)

        fill_queue_proc.join()


class PsiTsv(AnalysisTypeTsv):
    def __init__(self):
        super().__init__(ViewPsi)
        self.tab_output()

    def tsv_row(self, q, e, tsv_file, fieldnames):
        log = voila_log()
        with ViewPsi() as m, ViewSpliceGraph() as sg:

            with open(tsv_file, 'a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene_id in self.gene_ids(q, e):

                    gene = sg.gene(gene_id)

                    for lsv_id in m.view_gene_lsvs(gene.id):
                        lsv = m.lsv(lsv_id)
                        lsv_junctions = lsv.junctions
                        annot_juncs = sg.annotated_junctions(gene, lsv)
                        lsv_exons = sg.lsv_exons(gene, lsv)

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
                            'De Novo Junctions': semicolon_join(annot_juncs),
                            'Junctions coords': semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'Exons coords': semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in filter_exons(lsv_exons)
                            ),
                            # 'IR coords': self.semicolon_join(
                            #     '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                            # ),
                            'E(PSI) per LSV junction': semicolon_join(lsv.means),
                            'Var(E(PSI)) per LSV junction': semicolon_join(lsv.variances)
                        }

                        log.debug('Write TSV row for {0}'.format(lsv_id))

                        lock.acquire()
                        writer.writerow(row)
                        lock.release()

    def tab_output(self):

        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction',
                      'Var(E(PSI)) per LSV junction',
                      'LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions',
                      'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords']

        self.write_tsv(fieldnames)


class DeltaPsiTsv(AnalysisTypeTsv):
    def __init__(self):
        super().__init__(ViewDeltaPsi)
        self.tab_output()

    def tsv_row(self, q, e, tsv_file, fieldnames):
        log = voila_log()
        config = self.config

        with ViewDeltaPsi() as m, ViewSpliceGraph() as sg:
            metadata = m.view_metadata
            group1 = metadata['group_names'][0]
            group2 = metadata['group_names'][1]

            with open(tsv_file, 'a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene_id in self.gene_ids(q, e):

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
                            'De Novo Junctions': semicolon_join(annot_juncs),
                            'Junctions coords': semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'Exons coords': semicolon_join(
                                '{0}-{1}'.format(start, end) for start, end in filter_exons(lsv_exons)
                            ),
                            # TODO: this needs to be re-implemented
                            # 'IR coords': self.semicolon_join(
                            #     '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                            # ),
                            'E(dPSI) per LSV junction': semicolon_join(
                                excl_incl[i][1] - excl_incl[i][0] for i in
                                range(np.size(lsv.bins, 0))
                            ),
                            'P(|dPSI|>=%.2f) per LSV junction' % config.threshold: semicolon_join(
                                matrix_area(b, config.threshold) for b in lsv.bins
                            ),
                            'P(|dPSI|<=%.2f) per LSV junction' % config.non_changing_threshold: semicolon_join(
                                lsv.high_probability_non_changing()
                            ),
                            '%s E(PSI)' % group1: semicolon_join(
                                '%.3f' % i for i in group_means[group1]
                            ),
                            '%s E(PSI)' % group2: semicolon_join(
                                '%.3f' % i for i in group_means[group2]
                            )
                        }

                        lock.acquire()
                        writer.writerow(row)
                        lock.release()

    def tab_output(self):
        config = Config()

        with ViewDeltaPsi() as v:
            metadata = v.view_metadata

        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(dPSI) per LSV junction',
                      'P(|dPSI|>=%.2f) per LSV junction' % config.threshold,
                      'P(|dPSI|<=%.2f) per LSV junction' % config.non_changing_threshold,
                      '%s E(PSI)' % metadata['group_names'][0], '%s E(PSI)' % metadata['group_names'][1], 'LSV Type',
                      'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords']

        self.write_tsv(fieldnames)
