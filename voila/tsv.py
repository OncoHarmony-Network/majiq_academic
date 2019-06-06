import csv
import multiprocessing
import os
from datetime import datetime
from pathlib import Path
from queue import Empty

import numpy as np

from voila import constants
from voila.api.view_matrix import ViewHeterogens, ViewPsi, ViewDeltaPsi
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.config import ViewConfig, TsvConfig
from voila.exceptions import VoilaException, UnknownAnalysisType
from voila.view import views
from voila.vlsv import get_expected_psi, matrix_area
from voila.voila_log import voila_log

# lock used when writing files.
lock = multiprocessing.Lock()


def exon_str(exons):
    """
    Take in a list of exons and is an exon coordinate is -1, return 'nan' instead.  This helps with parsing the tsv file
    and something that is expected.
    :param exons: list of exon coordinates
    :return: exon coordinates or 'nan'
    """
    for exon in exons:
        if exon[0] == -1:
            yield 'nan', exon[1]
        elif exon[1] == -1:
            yield exon[0], 'nan'
        else:
            yield exon


def semicolon(value_list):
    """
    Take a list of value and return a semicolon separated list as a string.
    :param value_list: list of values
    :return: s=String
    """
    return ';'.join(str(x) for x in value_list)


def intron_retention_coords(lsv, juncs):
    """
    If this LSV has intron retention, then return IR coordinates with a hyphen separating them or return an empty
    string.  This is a helper function used in two classes.
    :param lsv:  LSV object
    :param juncs: list of LSV junctions which contain the coordinates of the IR at the end.
    :return: String
    """
    if lsv.alternative_intron:
        return '-'.join(map(str, juncs[-1]))
    else:
        return ''


class Tsv:
    def __init__(self):
        """
        Factory class used to create the TSV file for the specific analysis type.
        """
        config = TsvConfig()
        analysis_type = config.analysis_type

        voila_log().info(analysis_type + ' TSV')

        if analysis_type == constants.ANALYSIS_PSI:
            PsiTsv()
        elif analysis_type == constants.ANALYSIS_DELTAPSI:
            DeltaPsiTsv()
        elif analysis_type == constants.ANALYSIS_HETEROGEN:
            HeterogenTsv()
        else:
            raise UnknownAnalysisType(analysis_type)


class AnalysisTypeTsv:
    def __init__(self, view_matrix):
        """
        Super class to each of the analysis types classes for creating a TSV file.
        :param view_matrix: the view matrix class for the analysis type you're attempting to create a TSV file for.
        """

        start_time = datetime.now()
        config = TsvConfig()
        log = voila_log()

        self.filter_gene_ids = set()
        self.filter_lsv_ids = set()
        self.view_matrix = view_matrix

        if config.show_all:
            log.info('Showing all results and ignoring all filters. ')
        else:
            self.validate_filters()

        with view_matrix() as m:
            self.group_names = m.group_names

        self.tab_output()

        elapsed_time = datetime.now() - start_time
        voila_log().info('Duration: ' + str(elapsed_time))

    def tsv_row(self, q, e, tsv_file, fieldnames):
        """
        Used in multi-processing to write each row of the tsv file.
        :param q: Queue containing gene_ids
        :param e: Event flag for when the queue is no long having elements pushed to it.
        :param tsv_file: TSV filename.
        :param fieldnames: column fieldnames for the the TSV file.
        :return: None
        """

        raise NotImplementedError()

    def validate_filters(self):
        """
        Validate command line filters.  If filters are listed but not found, then they will printed out to the log.  If
        all filters are not found, then an error will be raised.
        :return: None
        """

        config = TsvConfig()
        log = voila_log()

        if config.lsv_ids:
            log.info('Validating LSV ids filter...')

            with self.view_matrix() as m:
                for lsv_id in config.lsv_ids:
                    if m.lsv(lsv_id).exists:
                        self.filter_lsv_ids.add(lsv_id)

            not_found_lsvs = set(config.lsv_ids) - self.filter_lsv_ids
            if not_found_lsvs:
                log.warning('LSV IDs not found: ' + ', '.join(not_found_lsvs))

        if config.lsv_types:
            log.info('Validating LSV types filter...')
            found_lsv_types = set()
            with self.view_matrix() as m:
                for lsv in m.lsvs():
                    lsv_type = lsv.lsv_type
                    if lsv_type in config.lsv_types:
                        self.filter_lsv_ids.add(lsv.lsv_id)
                        found_lsv_types.add(lsv_type)

            not_found_lsvs = set(config.lsv_types) - found_lsv_types
            if not_found_lsvs:
                log.warning('LSV types not found: ' + ', '.join(not_found_lsvs))

        if config.gene_ids:
            log.info('Validating gene IDs filter...')

            with self.view_matrix() as m:
                for gene_id in config.gene_ids:
                    if any(m.lsv_ids([gene_id])):
                        self.filter_gene_ids.add(gene_id)

            not_found_genes = set(config.gene_ids) - self.filter_gene_ids
            if not_found_genes:
                log.warning('Genes IDs not found: ' + ', '.join(not_found_genes))

        if config.gene_names:
            log.info('Validating gene names filter...')
            found_gene_names = set()
            with self.view_matrix() as m:
                matrix_gene_ids = list(m.gene_ids)

            with ViewSpliceGraph() as sg:
                for gene in sg.genes():
                    if gene['name'] in config.gene_names and gene['id'] in matrix_gene_ids:
                        self.filter_gene_ids.add(gene['id'])
                        found_gene_names.add(gene['name'])

            not_found_genes = set(config.gene_names) - found_gene_names
            if not_found_genes:
                log.warning('Gene names not found: ' + ', '.join(not_found_genes))

        if any((config.gene_names, config.gene_ids, config.lsv_ids, config.lsv_types)) and not any(
                (self.filter_gene_ids, self.filter_lsv_ids)):
            raise VoilaException('Filters were specified, but values were found in Voila or Splice Graph files.')

    def lsvs(self, gene_id):
        """
        Returns a generator that handles the filtering of LSVs based on the command line arguments.
        :param gene_id: Gene identifier for LSVs.
        :return: Generator
        """
        config = TsvConfig()

        with self.view_matrix() as m:
            lsvs = m.lsvs(gene_id)

            if config.show_all:

                yield from lsvs

            else:

                if self.filter_lsv_ids:
                    lsvs = list(lsv for lsv in lsvs if lsv.lsv_id in self.filter_lsv_ids)

                if self.filter_gene_ids:
                    lsvs = list(lsv for lsv in lsvs if lsv.gene_id in self.filter_gene_ids)

                if config.probability_threshold:
                    t = config.threshold
                    p = config.probability_threshold
                    lsvs = list(lsv for lsv in lsvs if any(matrix_area(b, t) >= p for b in lsv.bins))

                if config.analysis_type != constants.ANALYSIS_HETEROGEN:
                    lsvs = list(lsv for lsv in lsvs if any(m >= config.threshold for m in np.abs(lsv.means)))

                yield from lsvs

    def tab_output(self):
        """
        Helper method to set up fieldnames and then call self.write_tsv().
        :return: None
        """

        raise NotImplementedError()

    def fill_queue(self, queue, event):
        """
        Fill queue with gene IDs.
        :param queue: Queue to hold Gene IDs
        :param event: Event to flag when all Genes have been added to queue.
        :return: None
        """

        with self.view_matrix() as m:
            for gene_id in m.gene_ids:
                queue.put(gene_id)
            event.set()

    @staticmethod
    def gene_ids(q, e):
        """
        Get gene ids from queue.
        :param q: Queue containing gene IDS.
        :param e: Event flag to check if queue is still being added to.
        :return: None
        """

        while not (e.is_set() and q.empty()):
            try:
                yield q.get_nowait()
            except Empty:
                pass

        assert q.empty()

    def write_tsv(self, fieldnames):
        """
        Using the supplied fieldnames, start the pool writing the TSV file.
        :param fieldnames: list of column field names.
        :return: None
        """

        log = voila_log()
        log.info("Creating Tab-delimited output file")

        config = ViewConfig()
        nproc = config.nproc
        multiple_results = []

        tsv_file = config.file_name
        os.makedirs(os.path.dirname(tsv_file) or '.', exist_ok=True)
        tsv_file = Path(tsv_file)

        mgr = multiprocessing.Manager()

        queue = mgr.Queue(nproc * 2)
        event = mgr.Event()

        fill_queue_proc = multiprocessing.Process(target=self.fill_queue, args=(queue, event))
        fill_queue_proc.start()

        with tsv_file.open('w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

        with multiprocessing.Pool(processes=nproc) as pool:
            for _ in range(nproc):
                multiple_results.append(pool.apply_async(self.tsv_row, (queue, event, tsv_file, fieldnames)))
            [r.get() for r in multiple_results]

        fill_queue_proc.join()

        log.info("Delimited output file successfully created in: " + str(tsv_file))


class PsiTsv(AnalysisTypeTsv):
    def __init__(self):
        """
        Class to write TSV file for PSI analysis type.
        """

        super().__init__(ViewPsi)

    def tsv_row(self, q, e, tsv_file, fieldnames):
        log = voila_log()

        with ViewSpliceGraph() as sg:
            genome = sg.genome

            with tsv_file.open('a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene_id in self.gene_ids(q, e):

                    gene = sg.gene(gene_id)
                    chromosome = gene['chromosome']

                    for psi in self.lsvs(gene_id):
                        lsv_id = psi.lsv_id
                        lsv_junctions = psi.junctions
                        annot_juncs = sg.annotated_junctions(gene_id, lsv_junctions)
                        lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
                        ir_coords = intron_retention_coords(psi, lsv_junctions)
                        start, end = views.lsv_boundries(lsv_exons)

                        row = {
                            '#Gene Name': gene['name'],
                            'Gene ID': gene_id,
                            'LSV ID': lsv_id,
                            'LSV Type': psi.lsv_type,
                            'A5SS': psi.a5ss,
                            'A3SS': psi.a3ss,
                            'ES': psi.exon_skipping,
                            'Num. Junctions': psi.junction_count,
                            'Num. Exons': psi.exon_count,
                            'chr': gene['chromosome'],
                            'strand': gene['strand'],
                            'De Novo Junctions': semicolon(annot_juncs),
                            'Junctions coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'Exons coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in exon_str(lsv_exons)
                            ),
                            'IR coords': ir_coords,
                            'E(PSI) per LSV junction': semicolon(psi.means),
                            'Var(E(PSI)) per LSV junction': semicolon(psi.variances),
                            'UCSC LSV Link': views.ucsc_href(genome, chromosome, start, end)
                        }

                        lock.acquire()
                        log.debug('Write TSV row for {0}'.format(lsv_id))
                        writer.writerow(row)
                        lock.release()
                    q.task_done()

    def tab_output(self):
        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction',
                      'LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords', 'UCSC LSV Link']

        self.write_tsv(fieldnames)


class HeterogenTsv(AnalysisTypeTsv):
    def __init__(self):
        """
        Class to write TSV file for Heterogen analysis type.
        """

        super().__init__(ViewHeterogens)

    def tab_output(self):
        with ViewHeterogens() as m:
            group_names = m.group_names
            stats_column_names = list(m.junction_stats_column_names)

            fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'LSV Type', 'strand', 'chr'] + \
                         ['%s E(PSI)' % group for group in group_names] + stats_column_names + \
                         ['A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions',
                          'Junctions coords', 'Exons coords', 'IR coords', 'UCSC LSV Link']

        self.write_tsv(fieldnames)

    def tsv_row(self, q, e, tsv_file, fieldnames):
        log = voila_log()
        # config = TsvConfig()
        group_names = self.group_names

        with ViewSpliceGraph() as sg:
            genome = sg.genome

            with tsv_file.open('a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene_id in self.gene_ids(q, e):

                    gene = sg.gene(gene_id)
                    chromosome = gene['chromosome']

                    for het in self.lsvs(gene_id):
                        lsv_id = het.lsv_id

                        lsv_junctions = het.junctions
                        annot_juncs = sg.annotated_junctions(gene_id, lsv_junctions)
                        lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
                        mean_psi = list(het.mean_psi)
                        ir_coords = intron_retention_coords(het, lsv_junctions)
                        start, end = views.lsv_boundries(lsv_exons)

                        row = {
                            '#Gene Name': gene['name'],
                            'Gene ID': gene_id,
                            'LSV ID': lsv_id,
                            'LSV Type': het.lsv_type,
                            'A5SS': het.a5ss,
                            'A3SS': het.a3ss,
                            'ES': het.exon_skipping,
                            'Num. Junctions': len(lsv_junctions),
                            'Num. Exons': het.exon_count,
                            'chr': gene['chromosome'],
                            'strand': gene['strand'],
                            'De Novo Junctions': semicolon(annot_juncs),
                            'Junctions coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'Exons coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in exon_str(lsv_exons)
                            ),
                            'IR coords': ir_coords,
                            'UCSC LSV Link': views.ucsc_href(genome, chromosome, start, end)
                        }

                        for grp, mean in zip(group_names, np.array(mean_psi).transpose((1, 0, 2))):
                            row[grp + ' E(PSI)'] = semicolon(get_expected_psi(x) for x in mean)

                        for key, value in het.junction_stats:
                            row[key] = semicolon(value)

                        lock.acquire()
                        log.debug('Write TSV row for {0}'.format(lsv_id))
                        writer.writerow(row)
                        lock.release()

                    q.task_done()


class DeltaPsiTsv(AnalysisTypeTsv):
    def __init__(self):
        """
        Class to write TSV file for Delta PSI analysis type.
        """
        super().__init__(ViewDeltaPsi)

    def tsv_row(self, q, e, tsv_file, fieldnames):
        log = voila_log()
        config = TsvConfig()
        group1, group2 = self.group_names

        with ViewSpliceGraph() as sg:
            genome = sg.genome

            with tsv_file.open('a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                for gene_id in self.gene_ids(q, e):

                    gene = sg.gene(gene_id)
                    chromosome = gene['chromosome']

                    for dpsi in self.lsvs(gene_id):
                        lsv_id = dpsi.lsv_id

                        lsv_junctions = dpsi.junctions
                        annot_juncs = sg.annotated_junctions(gene_id, lsv_junctions)
                        lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
                        excl_incl = dpsi.excl_incl
                        group_means = dict(dpsi.group_means)
                        bins = dpsi.bins
                        ir_coords = intron_retention_coords(dpsi, lsv_junctions)
                        start, end = views.lsv_boundries(lsv_exons)

                        row = {
                            '#Gene Name': gene['name'],
                            'Gene ID': gene_id,
                            'LSV ID': lsv_id,
                            'LSV Type': dpsi.lsv_type,
                            'A5SS': dpsi.a5ss,
                            'A3SS': dpsi.a3ss,
                            'ES': dpsi.exon_skipping,
                            'Num. Junctions': dpsi.junction_count,
                            'Num. Exons': dpsi.exon_count,
                            'chr': gene['chromosome'],
                            'strand': gene['strand'],
                            'De Novo Junctions': semicolon(annot_juncs),
                            'Junctions coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'Exons coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in exon_str(lsv_exons)
                            ),
                            'IR coords': ir_coords,
                            'E(dPSI) per LSV junction': semicolon(
                                excl_incl[i][1] - excl_incl[i][0] for i in
                                range(np.size(bins, 0))
                            ),
                            'P(|dPSI|>=%.2f) per LSV junction' % config.threshold: semicolon(
                                matrix_area(b, config.threshold) for b in bins
                            ),
                            'P(|dPSI|<=%.2f) per LSV junction' % config.non_changing_threshold: semicolon(
                                dpsi.high_probability_non_changing()
                            ),
                            '%s E(PSI)' % group1: semicolon(
                                '%.3f' % i for i in group_means[group1]
                            ),
                            '%s E(PSI)' % group2: semicolon(
                                '%.3f' % i for i in group_means[group2]
                            ),
                            'UCSC LSV Link': views.ucsc_href(genome, chromosome, start, end)
                        }

                        lock.acquire()
                        log.debug('Write TSV row for {0}'.format(lsv_id))
                        writer.writerow(row)
                        lock.release()

                    q.task_done()

    def tab_output(self):
        config = ViewConfig()

        with ViewDeltaPsi() as v:
            grp_names = v.group_names

            fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(dPSI) per LSV junction',
                          'P(|dPSI|>=%.2f) per LSV junction' % config.threshold,
                          'P(|dPSI|<=%.2f) per LSV junction' % config.non_changing_threshold,
                          '%s E(PSI)' % grp_names[0], '%s E(PSI)' % grp_names[1], 'LSV Type', 'A5SS', 'A3SS', 'ES',
                          'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr', 'strand', 'Junctions coords',
                          'Exons coords', 'IR coords', 'UCSC LSV Link']

        self.write_tsv(fieldnames)
