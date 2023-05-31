import json
from itertools import chain


import numpy as np

from rna_voila import constants
from rna_voila.api.matrix_hdf5 import MatrixHdf5
from rna_voila.api.view_matrix import ViewHeterogens, ViewDeltaPsi, ViewPsi, ViewPsis, ViewHeterogen
from rna_voila.api.view_splice_graph import ViewSpliceGraph
from rna_voila.config import ViewConfig
from rna_voila.exceptions import UnknownAnalysisType, IndexNotFound, UnknownIndexFieldType, GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile
from rna_voila.vlsv import matrix_area
from rna_voila.voila_log import voila_log
from multiprocessing import Pool, Manager
import time
import hashlib
import os

lsv_filters = ['a5ss', 'a3ss', 'exon_skipping', 'target', 'source', 'binary', 'complex', 'intron_retention']
psi_keys = ['lsv_id', 'gene_id', 'gene_name'] + lsv_filters
dpsi_keys = ['lsv_id', 'gene_id', 'gene_name', 'excl_incl', 'dpsi_threshold', 'confidence_threshold'] + lsv_filters
het_keys = ['lsv_id', 'gene_id', 'gene_name', 'dpsi_threshold', 'stat_threshold'] + lsv_filters

skip_strict_indexing = False

class Index:
    def __init__(self, force_create=False, voila_files=None):
        """
        Factory class to generate the index for the supplied analysis type.
        """

        self.config = ViewConfig()
        analysis_type = self.config.analysis_type
        self.force_index = force_create or self.config.force_index
        self.index_voila_files = voila_files

        global skip_strict_indexing
        skip_strict_indexing = self.config.skip_type_indexing

        # The case where there's no analysis type, we're assuming this is splice graph only.
        if analysis_type:

            if analysis_type == constants.ANALYSIS_PSI:
                index = self._psi
            elif analysis_type == constants.ANALYSIS_DELTAPSI:
                index = self._deltapsi
            elif analysis_type == constants.ANALYSIS_HETEROGEN:
                index = self._heterogen
            else:
                raise UnknownAnalysisType(analysis_type)

            index()

    @staticmethod
    def _index_in_voila(voila_file, remove_index=False):
        """
        Check if index has already been created in voila file. If remove_index has been set, then attempt to use the
        h5py api to remove the index dataset from the file.

        If the case of multiple voila files, we hash all of them and check that the hash of that group of hashes
        matches what is stored in the voila file. Here we are assuming that the 'order' of the files is not important
        (for example, if we check for a match in the first file but the previous index was stored in a different
        file, the index will be rebuilt again for the same data)
        :param voila_file:
        :param remove_index:
        :return:
        """

        if os.path.exists(voila_file):
            with MatrixHdf5(voila_file, 'r', pre_config=True) as m:
                return m.has_index()
        return False
        #
        # with h5py.File(voila_file, 'a') as h:
        #     index_in_h = 'index' in h
        #
        #     if remove_index and index_in_h:
        #         voila_log().info('Removing index from HDF5')
        #         del h['index']
        #
        #     # voila_files = ViewConfig().voila_files
        #
        #     # if 'input_hash' in h:
        #     #     prior = h.get('input_hash')[0].decode('utf-8')
        #     #     new = Index._get_files_hash(voila_files)
        #     #     index_in_h = (prior == new)
        #
        #     return index_in_h

    @staticmethod
    def _get_files_hash(voila_files):
        if ViewConfig().index_file:
            # if we use a separate file for indexing, we can verify the hash of all inputs in the verification
            h = hashlib.sha1()
            for filename in sorted(voila_files):
                with open(filename, 'rb') as f:
                    h.update(f.read())
            return h.hexdigest()
        else:
            # for now we get a hash based on the combined name of all files
            # (this is because we can not easily do it based on the content because we change the content for this process)
            return hashlib.sha1(''.join(voila_files).encode('utf-8')).hexdigest()

    @staticmethod
    def _write_index(voila_file, voila_index, dtype):
        """
        Helper method to write voila index to voila file using specific numpy data type.

        :param voila_file: location and name of voila file.
        :param voila_index: Array of index data.
        :param dtype: numpy data type string.
        :return: None
        """
        voila_index = np.array(voila_index, dtype=np.dtype(dtype))
        voila_files = ViewConfig().voila_files

        with MatrixHdf5(voila_file, 'a', pre_config=True) as m:
            hashval = Index._get_files_hash(voila_files)
            m.write_index(voila_index, hashval)


    @staticmethod
    def _get_voila_index_file():
        c = ViewConfig()
        if c.index_file:
            return c.index_file
        else:
            return c.voila_file

    @staticmethod
    def _create_dtype(voila_index):
        first_row = voila_index[0]
        dtype = []

        for field in first_row:
            if isinstance(field, str):
                dtype.append(len(field))
            elif isinstance(field, bool):
                dtype.append('?')
            elif isinstance(field, np.float64):
                dtype.append('f4')
            else:
                raise UnknownIndexFieldType(field)

        for row in voila_index:
            for idx, field in enumerate(row):
                if dtype[idx] not in ('?', 'f4'):
                    dtype[idx] = max(len(field), dtype[idx])

        dtype = [d if d in ('?', 'f4') else 'S' + str(d) for d in dtype]
        dtype = ','.join(dtype)

        return dtype

    @staticmethod
    def _heterogen_pool_add_index(args):
        """
        Multithread inner function for each iteration of _heterogen loop below
        """
        lsv_id, g, vfs, q = args
        config = ViewConfig()

        g_dpsi_thresh = np.array([-1])

        with ViewHeterogens() as m:
            het = m.lsv(lsv_id)
            gene_id = het.gene_id
            gene_name = g[gene_id]
            g_stats_thresh = np.array([2 for _ in m.stat_names])

            if skip_strict_indexing:
                lsv_f = [True for _ in lsv_filters]
            else:
                lsv_f = [getattr(het, f) for f in lsv_filters]

        for vf in config.voila_files:
            if vfs and vf not in vfs:
                continue
            with ViewHeterogen(vf) as m:
                try:
                    het = m.lsv(lsv_id)

                    dpsi_thresh = het.dpsi_median()
                    g_dpsi_thresh = np.maximum(np.max(np.abs(dpsi_thresh)), g_dpsi_thresh)

                    # get min val for each stat
                    stats_thresh = het.junction_stats
                    g_stats_thresh = np.minimum(np.min(stats_thresh, axis=0), g_stats_thresh)
                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    continue

        g_dpsi_thresh = g_dpsi_thresh.tolist()[0]
        g_dpsi_thresh = json.dumps(g_dpsi_thresh)

        g_stats_thresh = g_stats_thresh.tolist()
        g_stats_thresh = json.dumps(g_stats_thresh)

        row = (lsv_id, gene_id, gene_name, g_dpsi_thresh, g_stats_thresh)

        # For some reason, numpy needs these in tuples.
        row = tuple(chain(row, lsv_f))

        if q:
            q.put(row)
        return row

    def _heterogen(self):
        """
        Create index for heterogen analysis type.  This is an odd case as there can be more then one voila file.  We
        create the index in the first voila file.  First is determined by ordering the file name and location.

        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = self.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            with ViewSpliceGraph() as sg:
                g = sg.gene_ids2gene_names

            if config.nproc > 1:
                manager = Manager()
                q = manager.Queue()

                with ViewHeterogens() as m:
                    lsv_ids = [(x, g, self.index_voila_files, q) for x in m.lsv_ids()]


                p = Pool(config.nproc)
                work_size = len(lsv_ids)

                # voila_index = p.map(self._heterogen_pool_add_index, zip(lsv_ids, range(work_size), repeat(work_size)))
                voila_index = p.map_async(self._heterogen_pool_add_index, lsv_ids)

                # monitor loop
                while True:
                    if voila_index.ready():
                        break
                    else:
                        size = q.qsize()
                        print("Indexing LSV IDs: %d / %d" % (size, work_size))
                        time.sleep(2)

                log.info('Writing index: ' + voila_file)
                voila_index = voila_index.get()

            else:
                voila_index = []
                with ViewHeterogens() as m:
                    for x in m.lsv_ids():
                        args = (x, g, self.index_voila_files, None)
                        voila_index.append(self._heterogen_pool_add_index(args))

            dtype = self._create_dtype(voila_index)
            self._write_index(voila_file, voila_index, dtype)
        else:
            log.info('Using index: ' + voila_file)

    @staticmethod
    def _deltapsi_pool_add_index(args):
        """
        Multithread inner function for each iteration of _deltapsi loop below
        """
        lsv_id, g, q = args

        with ViewDeltaPsi() as m:
            dpsi = m.lsv(lsv_id)
            excl_incl = max(abs(a - b) for a, b in dpsi.excl_incl)
            gene_id = dpsi.gene_id

            dpsi_thresh = dpsi.means
            dpsi_thresh = np.abs(dpsi_thresh)
            dpsi_thresh = dpsi_thresh.tolist()
            dpsi_thresh = json.dumps(dpsi_thresh)

            # Calculate a list of confidence for 10 (0 thru 1) values of threshold. This is done because we
            # don't want to store the bins data in the index.  Although, this might be a good idea once the
            # index gets large enough.
            bins = dpsi.bins
            confidence_thresh = list(max(matrix_area(b, x) for b in bins) for x in np.linspace(0, 1, 10))
            confidence_thresh = json.dumps(confidence_thresh)

            gene_name = g[gene_id]

            row = (lsv_id, gene_id, gene_name, excl_incl, dpsi_thresh, confidence_thresh)

            if skip_strict_indexing:
                lsv_f = [True for _ in lsv_filters]
            else:
                lsv_f = [getattr(dpsi, f) for f in lsv_filters]


        # For some reason, numpy needs these in tuples.
        row = tuple(chain(row, lsv_f))
        if q:
            q.put(row)
        return row

    def _deltapsi(self):
        """
        Generates index for delta psi analysis type.

        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = self.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            with ViewSpliceGraph() as sg:
                g = sg.gene_ids2gene_names

            if config.nproc > 1:

                manager = Manager()
                q = manager.Queue()

                with ViewDeltaPsi() as m:
                    lsv_ids = [(x, g, q) for x in m.lsv_ids()]
                p = Pool(config.nproc)
                work_size = len(lsv_ids)

                voila_index = p.map_async(self._deltapsi_pool_add_index, lsv_ids)

                # monitor loop
                while True:
                    if voila_index.ready():
                        break
                    else:
                        size = q.qsize()
                        print("Indexing LSV IDs: %d / %d" % (size, work_size))
                        time.sleep(2)

                log.info('Writing index: ' + voila_file)
                voila_index = voila_index.get()

            else:
                voila_index = []
                with ViewDeltaPsi() as m:
                    for x in m.lsv_ids():
                        args = (x, g, None)
                        voila_index.append(self._deltapsi_pool_add_index(args))

            dtype = self._create_dtype(voila_index)
            self._write_index(voila_file, voila_index, dtype)
        else:
            log.info('Using index: ' + voila_file)

    @staticmethod
    def _psi_pool_add_index(args):
        """
        Multithread inner function for each iteration of _psi loop below
        """
        lsv_id, g, m, q = args
        lsv = m.lsv(lsv_id)
        gene_id = lsv.gene_id

        gene_name = g[gene_id]

        row = (lsv_id, gene_id, gene_name)

        if skip_strict_indexing:
            lsv_f = [True for _ in lsv_filters]
        else:
            lsv_f = [getattr(lsv, f) for f in lsv_filters]

        # For some reason, numpy needs these in tuples.
        row = tuple(chain(row, lsv_f))

        if q:
            q.put(row)
        return row

    def _psi(self):
        """
        Create index to PSI analysis type.
        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = self.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            with ViewSpliceGraph() as sg:
                g = sg.gene_ids2gene_names

            if config.nproc > 1:

                manager = Manager()
                q = manager.Queue()

                with ViewPsis() as m:
                    lsv_ids = [(x, g, m, q) for x in m.lsv_ids()]
                p = Pool(config.nproc)
                work_size = len(lsv_ids)

                voila_index = p.map_async(self._psi_pool_add_index, lsv_ids)

                # monitor loop
                while True:
                    if voila_index.ready():
                        break
                    else:
                        size = q.qsize()
                        print("Indexing LSV IDs: %d / %d" % (size, work_size))
                        time.sleep(2)

                log.info('Writing index: ' + voila_file)
                voila_index = voila_index.get()

            else:
                voila_index = []
                with ViewPsis() as m:
                    for x in m.lsv_ids():
                        args = (x, g, m, None)
                        voila_index.append(self._psi_pool_add_index(args))

            dtype = self._create_dtype(voila_index)
            self._write_index(voila_file, voila_index, dtype)
        else:
            log.info('Using index: ' + voila_file)

    @staticmethod
    def _row_data(gene_id, keys):
        """
        For each row in index, zip list of keys with values in the row.
        :param gene_id: gene id
        :param keys: index field names
        :return:
        """

        index_file = Index._get_voila_index_file()

        try:
            gene_id = gene_id.encode('utf-8')
        except AttributeError:
            pass

        with MatrixHdf5(index_file, 'r') as m:

            try:
                for row in m.get_index():
                    if gene_id is None or gene_id == row[1]:
                        yield dict(zip(keys, row))

            except KeyError:
                raise IndexNotFound()

    @classmethod
    def psi(cls, gene_id=None):
        """
        Get PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        yield from cls._row_data(gene_id, psi_keys)

    @classmethod
    def delta_psi(cls, gene_id=None):
        """
        Get Delta PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        yield from cls._row_data(gene_id, dpsi_keys)

    @classmethod
    def heterogen(cls, gene_id=None):
        """
        Get Heterogen index data in a dictionary for each row.
        :return:
        """

        yield from cls._row_data(gene_id, het_keys)
