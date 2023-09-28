import csv
import multiprocessing
import os, sys
from datetime import datetime
from pathlib import Path
from queue import Empty

import numpy as np
import json

from rna_voila import constants
from rna_voila.api.view_matrix import ViewHeterogens, ViewPsi, ViewDeltaPsi
from rna_voila.api.view_splice_graph import ViewSpliceGraph
from rna_voila.config import TsvConfig
from rna_voila.exceptions import VoilaException, UnknownAnalysisType
from rna_voila.view import views
from rna_voila.vlsv import get_expected_psi, matrix_area
from rna_voila.voila_log import voila_log
from rna_voila.api.view_matrix import ViewMatrix

# lock used when writing files.
lock = multiprocessing.Lock()


def exon_str(exons):
    """
    Take in a list of exons and is an exon coordinate is -1, return 'na' instead.  This helps with parsing the tsv file
    and something that is expected.
    :param exons: list of exon coordinates
    :return: exon coordinates or 'na'
    """
    for exon in exons:
        if exon[0] == -1:
            yield 'na', exon[1]
        elif exon[1] == -1:
            yield exon[0], 'na'
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
    if lsv.intron_retention:
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

        m_all = ViewMatrix()
        warnings = m_all.check_group_consistency()
        if warnings:
            for warning in warnings:
                voila_log().warning(f'Warning: detected groups with the same name "{warning[0]}", which have different sets of experiments: {warning[1]}')
            if not config.ignore_inconsistent_group_errors:
                voila_log().critical("Exiting due to previous warnings, pass --ignore-inconsistent-group-errors to run anyway")
                sys.exit(1)

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
            self._experiment_names = m.experiment_names
            self.experiment_names = []
            for group in self._experiment_names:
                for expname in group:
                    if expname:
                        self.experiment_names.append(expname)

        self.tab_output()

        elapsed_time = datetime.now() - start_time
        voila_log().info('Duration: ' + str(elapsed_time))

    def get_metadata(self):
        return self.get_base_metadata()

    def get_base_metadata(self):
        with self.view_matrix() as m:
            metadata = {
                'voila_version': constants.VERSION,
                'command': ' '.join(sys.argv),
                'group_sizes': {
                    grp: sum(bool(x) for x in names)  # skip empty names used for padding
                    for grp, names in zip(m.group_names, m.experiment_names)
                },
            }

        return metadata

    def tsv_row(self, q, e, tsv_file, fieldnames, gene_ids=None):
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

        config = TsvConfig()
        nproc = config.nproc
        multiple_results = []

        tsv_file = config.file_name
        os.makedirs(os.path.dirname(tsv_file) or '.', exist_ok=True)
        tsv_file = Path(tsv_file)

        with tsv_file.open('w') as tsv:

            metadata = self.get_metadata()

            metadata_string = '\n'.join([f"# {l}" for l in json.dumps(metadata, sort_keys=True,
                                        indent=4, separators=(',', ': ')).split('\n')]) + '\n'
            tsv.write(metadata_string)

            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

        if nproc > 1:
            mgr = multiprocessing.Manager()

            queue = mgr.Queue(nproc * 2)
            event = mgr.Event()

            fill_queue_proc = multiprocessing.Process(target=self.fill_queue, args=(queue, event))
            fill_queue_proc.start()

            with multiprocessing.Pool(processes=nproc) as pool:
                for _ in range(nproc):
                    multiple_results.append(pool.apply_async(self.tsv_row, (queue, event, tsv_file, fieldnames)))
                [r.get() for r in multiple_results]

            fill_queue_proc.join()
        else:
            with self.view_matrix() as m:
                self.tsv_row(None, None, tsv_file, fieldnames, m.gene_ids)



        log.info("Delimited output file successfully created in: " + str(tsv_file))


class PsiTsv(AnalysisTypeTsv):
    def __init__(self):
        """
        Class to write TSV file for PSI analysis type.
        """

        super().__init__(ViewPsi)

    def get_metadata(self):
        return self.get_base_metadata()

    def tsv_row(self, q, e, tsv_file, fieldnames, gene_ids=None):
        log = voila_log()

        with ViewSpliceGraph() as sg:
            genome = sg.genome

            with tsv_file.open('a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                _gene_ids = list(self.gene_ids(q, e) if q else gene_ids)
                work_size = len(_gene_ids)

                for i, gene_id in enumerate(_gene_ids):

                    gene = sg.gene(gene_id)
                    chromosome = gene['chromosome']

                    if i % 10 == 0:
                        print('Processing rows [%d/%d]\r' % (i, work_size), end="")

                    for psi in self.lsvs(gene_id):
                        lsv_id = psi.lsv_id
                        lsv_junctions = psi.junctions
                        annot_juncs = sg.annotated_junctions(gene_id, lsv_junctions)
                        lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
                        ir_coords = intron_retention_coords(psi, lsv_junctions)
                        start, end = views.lsv_boundries(lsv_exons)



                        row = {
                            'gene_name': gene['name'],
                            'gene_id': gene_id,
                            'lsv_id': lsv_id,
                            'lsv_type': psi.lsv_type,
                            'num_junctions': psi.junction_count,
                            'num_exons': psi.exon_count,
                            'seqid': gene['chromosome'],
                            'strand': gene['strand'],
                            'de_novo_junctions': semicolon((1 if x == 0 else 0 for x in annot_juncs)),  # reverse flags
                            'junctions_coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'exons_coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in exon_str(lsv_exons)
                            ),
                            'ir_coords': ir_coords,
                            'mean_psi_per_lsv_junction': semicolon(f'{x:0.4f}' for x in psi.means),
                            'stdev_psi_per_lsv_junction': semicolon(f'{x:0.4f}' for x in psi.standard_deviations),
                            'ucsc_lsv_link': views.ucsc_href(genome, chromosome, start, end)
                        }

                        config = TsvConfig()
                        if config.show_read_counts:
                            experiment_reads = sg.lsv_reads(gene_id, lsv_junctions, self.experiment_names,
                                                            bool(ir_coords))
                            for exp in experiment_reads:
                                if exp in self.experiment_names:
                                    junc_reads, int_reads = experiment_reads[exp]
                                    row[f"{exp}_junction_reads"] = semicolon(junc_reads)
                                    row[f"{exp}_intron_retention_reads"] = semicolon(int_reads)

                        if lock:
                            lock.acquire()
                        log.debug('Write TSV row for {0}'.format(lsv_id))
                        writer.writerow(row)
                        if lock:
                            lock.release()
                    if q:
                        q.task_done()

        print('                                                  \r', end="")

    def tab_output(self):
        fieldnames = ['gene_name', 'gene_id', 'lsv_id', 'mean_psi_per_lsv_junction', 'stdev_psi_per_lsv_junction',
                      'lsv_type', 'num_junctions', 'num_exons', 'de_novo_junctions', 'seqid',
                      'strand', 'junctions_coords', 'exons_coords', 'ir_coords', 'ucsc_lsv_link']

        config = TsvConfig()
        if config.show_read_counts:
            for exp in self.experiment_names:
                fieldnames.append(f"{exp}_junction_reads")
                fieldnames.append(f"{exp}_intron_retention_reads")


        self.write_tsv(fieldnames)


class HeterogenTsv(AnalysisTypeTsv):
    def __init__(self):
        """
        Class to write TSV file for Heterogen analysis type.
        """

        self._quantiles = (0.25, 0.75,)
        super().__init__(ViewHeterogens)

    def get_metadata(self):
        metadata = self.get_base_metadata()
        with ViewHeterogens() as m:
            metadata['stat_names'] = m.stat_names
            metadata['psi_samples'] = m.psi_samples_summary
            metadata['test_percentile'] = m.test_percentile_summary
        return metadata

    def tab_output(self):
        with ViewHeterogens() as m:
            group_names = m.group_names
            stats_column_names = (
                list(m.junction_stats_column_names)
                + list(m.junction_psisamples_stats_column_names)
                + list(m.junction_scores_column_names)
            )

            fieldnames = [
                'gene_name',
                'gene_id',
                'lsv_id',
                'lsv_type',
                'strand',
                'seqid',
                *(
                    f'{group}_median_psi' for group in group_names
                ),
                *(
                    f'{group}_percentile{quant * 100:02.0f}_psi' for group in group_names for quant in self._quantiles
                ),
                *(
                    f'{group}_num_quantified' for group in group_names
                ),
                *stats_column_names,
                *m.changing_column_names,
                *m.nonchanging_column_names,
                'num_junctions',
                'num_exons',
                'de_novo_junctions',
                'junctions_coords',
                'exons_coords',
                'ir_coords',
                'ucsc_lsv_link'
            ]

        config = TsvConfig()
        if config.show_read_counts:
            for exp in self.experiment_names:
                fieldnames.append(f"{exp}_junction_reads")
                fieldnames.append(f"{exp}_intron_retention_reads")

        if config.show_per_sample_psi:
            fieldnames += [f'{exp}_psi' for i in range(len(group_names)) for exp in self._experiment_names[i] ]

        self.write_tsv(fieldnames)

    def tsv_row(self, q, e, tsv_file, fieldnames, gene_ids=None):
        log = voila_log()
        config = TsvConfig()
        group_names = self.group_names

        with ViewSpliceGraph() as sg:
            genome = sg.genome

            with tsv_file.open('a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                _gene_ids = list(self.gene_ids(q, e) if q else gene_ids)
                work_size = len(_gene_ids)

                for i, gene_id in enumerate(_gene_ids):

                    gene = sg.gene(gene_id)
                    chromosome = gene['chromosome']

                    if i % 10 == 0:
                        print('Processing rows [%d/%d]\r' % (i, work_size), end="")

                    for het in self.lsvs(gene_id):
                        lsv_id = het.lsv_id

                        lsv_junctions = het.junctions
                        annot_juncs = sg.annotated_junctions(gene_id, lsv_junctions)
                        lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
                        ir_coords = intron_retention_coords(het, lsv_junctions)
                        start, end = views.lsv_boundries(lsv_exons)

                        row = {
                            'gene_name': gene['name'],
                            'gene_id': gene_id,
                            'lsv_id': lsv_id,
                            'lsv_type': het.lsv_type,
                            'num_junctions': len(lsv_junctions),
                            'num_exons': het.exon_count,
                            'seqid': gene['chromosome'],
                            'strand': gene['strand'],
                            'de_novo_junctions': semicolon((1 if x == 0 else 0 for x in annot_juncs)),  # reverse flags
                            'junctions_coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'exons_coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in exon_str(lsv_exons)
                            ),
                            'ir_coords': ir_coords,
                            'ucsc_lsv_link': views.ucsc_href(genome, chromosome, start, end),
                            **{key: semicolon(values) for key, values in
                               het.changing(config.changing_pvalue_threshold,
                                            config.changing_between_group_dpsi)},
                            **{key: semicolon(values) for key, values in
                               het.nonchanging(config.non_changing_pvalue_threshold,
                                               config.non_changing_within_group_iqr,
                                               config.non_changing_between_group_dpsi)},
                        }

                        for grp, medians in zip(group_names, het.median_psi):
                            if (medians < 0).all():
                                row[f'{grp}_median_psi'] = 'NA'
                            else:
                                row[f'{grp}_median_psi'] = semicolon(f'{x:0.4f}' for x in medians)

                        for grp, num_quantified in zip(group_names, het.groups_quantified):
                            row[f'{grp}_num_quantified'] = num_quantified

                        for quant in self._quantiles:
                            for grp, medians in zip(group_names, het.quantile_psi(quant)):
                                if (medians < 0).all():
                                    row[f'{grp}_percentile{quant * 100:02.0f}_psi'] = 'NA'
                                else:
                                    row[f'{grp}_percentile{quant * 100:02.0f}_psi'] = semicolon(f'{x:0.3e}' for x in medians)

                        for key, values in het.junction_stats:
                            row[key] = semicolon(f'{x:0.3e}' for x in values)

                        for key, values in het.junction_psisamples_stats:
                            row[key] = semicolon(f'{x:0.3e}' for x in values)

                        for key, values in het.junction_scores:
                            row[key] = semicolon(f'{x:0.3e}' for x in values)

                        config = TsvConfig()
                        if config.show_read_counts:
                            experiment_reads = sg.lsv_reads(gene_id, lsv_junctions, self.experiment_names,
                                                            bool(ir_coords))
                            for exp in experiment_reads:
                                if exp in self.experiment_names:
                                    junc_reads, int_reads = experiment_reads[exp]
                                    row[f"{exp}_junction_reads"] = semicolon(junc_reads)
                                    row[f"{exp}_intron_retention_reads"] = semicolon(int_reads)

                        if config.show_per_sample_psi:
                            mu_psi = het.mu_psi
                            for i, grp in enumerate(group_names):
                                for j, exp in enumerate(self._experiment_names[i]):
                                    try:
                                        row[f"{exp}_psi"] = semicolon(mu_psi[x][i][j] for x in range(len(lsv_junctions)))
                                    except IndexError:
                                        row[f"{exp}_psi"] = ''



                        if lock:
                            lock.acquire()
                        log.debug('Write TSV row for {0}'.format(lsv_id))
                        writer.writerow(row)
                        if lock:
                            lock.release()
                    if q:
                        q.task_done()

        print('                                                  \r', end="")

class DeltaPsiTsv(AnalysisTypeTsv):
    def __init__(self):
        """
        Class to write TSV file for Delta PSI analysis type.
        """
        super().__init__(ViewDeltaPsi)

    def get_metadata(self):
        metadata = self.get_base_metadata()
        config = TsvConfig()
        metadata['changing_threshold'] = config.threshold
        metadata['non_changing_threshold'] = config.non_changing_between_group_dpsi
        return metadata


    def tsv_row(self, q, e, tsv_file, fieldnames, gene_ids=None):
        log = voila_log()
        config = TsvConfig()
        group1, group2 = self.group_names

        with ViewSpliceGraph() as sg:
            genome = sg.genome

            with tsv_file.open('a') as tsv:
                writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')

                _gene_ids = list(self.gene_ids(q, e) if q else gene_ids)
                work_size = len(_gene_ids)

                for i, gene_id in enumerate(_gene_ids):

                    gene = sg.gene(gene_id)
                    chromosome = gene['chromosome']

                    if i % 10 == 0:
                        print('Processing rows [%d/%d]\r' % (i, work_size), end="")

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
                            'gene_name': gene['name'],
                            'gene_id': gene_id,
                            'lsv_id': lsv_id,
                            'lsv_type': dpsi.lsv_type,
                            'num_junctions': dpsi.junction_count,
                            'num_exons': dpsi.exon_count,
                            'seqid': gene['chromosome'],
                            'strand': gene['strand'],
                            'de_novo_junctions': semicolon((1 if x == 0 else 0 for x in annot_juncs)),  # reverse flags
                            'junctions_coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in lsv_junctions
                            ),
                            'exons_coords': semicolon(
                                '{0}-{1}'.format(start, end) for start, end in exon_str(lsv_exons)
                            ),
                            'ir_coords': ir_coords,
                            'mean_dpsi_per_lsv_junction': semicolon(
                                f'{excl_incl[i][1] - excl_incl[i][0]:0.4f}' for i in
                                range(np.size(bins, 0))
                            ),
                            'probability_changing': semicolon(
                                f'{matrix_area(b, config.threshold):0.3e}' for b in bins
                            ),
                            'probability_non_changing': semicolon(
                                f'{x:0.3e}' for x in dpsi.high_probability_non_changing()
                            ),
                            '%s_mean_psi' % group1: semicolon(
                                f'{x:0.4f}' for x in group_means[group1]
                            ),
                            '%s_mean_psi' % group2: semicolon(
                                f'{x:0.4f}' for x in group_means[group2]
                            ),
                            'ucsc_lsv_link': views.ucsc_href(genome, chromosome, start, end)
                        }

                        config = TsvConfig()
                        if config.show_read_counts:
                            experiment_reads = sg.lsv_reads(gene_id, lsv_junctions, self.experiment_names,
                                                            bool(ir_coords))
                            for exp in experiment_reads:
                                if exp in self.experiment_names:
                                    junc_reads, int_reads = experiment_reads[exp]
                                    row[f"{exp}_junction_reads"] = semicolon(junc_reads)
                                    row[f"{exp}_intron_retention_reads"] = semicolon(int_reads)

                        if lock:
                            lock.acquire()
                        log.debug('Write TSV row for {0}'.format(lsv_id))
                        writer.writerow(row)
                        if lock:
                            lock.release()
                    if q:
                        q.task_done()

        print('                                                  \r', end="")

    def tab_output(self):

        with ViewDeltaPsi() as v:
            grp_names = v.group_names

            fieldnames = ['gene_name', 'gene_id', 'lsv_id', 'mean_dpsi_per_lsv_junction',
                          'probability_changing',
                          'probability_non_changing',
                          '%s_mean_psi' % grp_names[0], '%s_mean_psi' % grp_names[1], 'lsv_type',
                          'num_junctions', 'num_exons', 'de_novo_junctions', 'seqid', 'strand', 'junctions_coords',
                          'exons_coords', 'ir_coords', 'ucsc_lsv_link']

        config = TsvConfig()
        if config.show_read_counts:
            for exp in self.experiment_names:
                fieldnames.append(f"{exp}_junction_reads")
                fieldnames.append(f"{exp}_intron_retention_reads")

        self.write_tsv(fieldnames)
