import json
from itertools import chain

import h5py
import numpy as np

from voila import constants
from voila.api.view_matrix import ViewHeterogens, ViewDeltaPsi, ViewPsi, ViewPsis
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.config import ViewConfig
from voila.exceptions import UnknownAnalysisType, IndexNotFound, UnknownIndexFieldType
from voila.vlsv import matrix_area
from voila.voila_log import voila_log
from multiprocessing import Pool, Manager
import time
import hashlib

lsv_filters = ['a5ss', 'a3ss', 'exon_skipping', 'target', 'source', 'binary', 'complex']
psi_keys = ['lsv_id', 'gene_id', 'gene_name'] + lsv_filters
dpsi_keys = ['lsv_id', 'gene_id', 'gene_name', 'excl_incl', 'dpsi_threshold', 'confidence_threshold'] + lsv_filters
het_keys = ['lsv_id', 'gene_id', 'gene_name'] + lsv_filters


class Index:
    def __init__(self):
        """
        Factory class to generate the index for the supplied analysis type.
        """

        analysis_type = ViewConfig().analysis_type

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

        with h5py.File(voila_file, 'a') as h:
            index_in_h = 'index' in h

            if remove_index and index_in_h:
                voila_log().info('Removing index from HDF5')
                del h['index']

            voila_files = ViewConfig().voila_files

            if 'input_hash' in h:
                prior = h.get('input_hash')[0].decode('utf-8')
                new = Index._get_files_hash(voila_files)
                index_in_h = (prior == new)

            return index_in_h

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

        with h5py.File(voila_file, 'a') as h:
            if 'index' in h:
                del h['index']
            if 'input_hash' in h:
                del h['input_hash']
            h.create_dataset('index', voila_index.shape, data=voila_index)
            hashval = Index._get_files_hash(voila_files)
            h.create_dataset("input_hash", (1,), dtype="S40", data=(hashval.encode('utf-8'),))

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
        with ViewSpliceGraph() as sg, ViewHeterogens() as m:
            lsv_id, q = args
            het = m.lsv(lsv_id)

            gene_id = het.gene_id
            gene = sg.gene(gene_id)
            gene_name = gene['name']

            row = (lsv_id, gene_id, gene_name)

            lsv_f = [getattr(het, f) for f in lsv_filters]

            # For some reason, numpy needs these in tuples.
            row = tuple(chain(row, lsv_f))
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
        force_index = remove_index = config.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            m = Manager()
            q = m.Queue()

            with ViewHeterogens() as m:
                lsv_ids = [(x, q) for x in m.lsv_ids()]
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

            dtype = self._create_dtype(voila_index)
            self._write_index(voila_file, voila_index, dtype)
        else:
            log.info('Using index: ' + voila_file)

    @staticmethod
    def _deltapsi_pool_add_index(args):
        """
        Multithread inner function for each iteration of _deltapsi loop below
        """
        with ViewSpliceGraph() as sg, ViewDeltaPsi() as m:
            lsv_id, q = args
            dpsi = m.lsv(lsv_id)

            gene_id = dpsi.gene_id
            gene = sg.gene(gene_id)
            gene_name = gene['name']

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

            excl_incl = dpsi.excl_incl
            excl_incl = max(abs(a - b) for a, b in excl_incl)

            lsv_f = [getattr(dpsi, f) for f in lsv_filters]
            row = (lsv_id, gene_id, gene_name, excl_incl, dpsi_thresh, confidence_thresh)

            # For some reason, numpy needs these in tuples.
            row = tuple(chain(row, lsv_f))
            q.put(row)
            return row

    def _deltapsi(self):
        """
        Generates index for delta psi analysis type.

        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = config.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            m = Manager()
            q = m.Queue()

            with ViewDeltaPsi() as m:
                lsv_ids = [(x, q) for x in m.lsv_ids()]
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

            dtype = self._create_dtype(voila_index)
            self._write_index(voila_file, voila_index, dtype)
        else:
            log.info('Using index: ' + voila_file)

    @staticmethod
    def _psi_pool_add_index(args):
        """
        Multithread inner function for each iteration of _psi loop below
        """
        with ViewSpliceGraph() as sg, ViewPsis() as m:
            lsv_id, q = args
            lsv = m.lsv(lsv_id)

            gene_id = lsv.gene_id
            gene = sg.gene(gene_id)
            gene_name = gene['name']

            row = (lsv_id, gene_id, gene_name)

            lsv_f = [getattr(lsv, f) for f in lsv_filters]

            # For some reason, numpy needs these in tuples.
            row = tuple(chain(row, lsv_f))
            q.put(row)
            return row

    def _psi(self):
        """
        Create index to PSI analysis type.
        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = config.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            m = Manager()
            q = m.Queue()

            with ViewPsis() as m:
                lsv_ids = [(x, q) for x in m.lsv_ids()]
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

        with h5py.File(index_file, 'r') as h:

            try:

                for row in h['index'].value:
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
