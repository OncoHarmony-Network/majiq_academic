import json
import os
from itertools import chain
from pathlib import Path

import h5py
import numpy as np

from voila import constants
from voila.api.view_matrix import ViewHeterogens, ViewDeltaPsi, ViewPsi
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.config import ViewConfig
from voila.exceptions import UnknownAnalysisType, IndexNotFound, UnknownIndexFieldType
from voila.vlsv import matrix_area
from voila.voila_log import voila_log

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

        :param voila_file:
        :param remove_index:
        :return:
        """

        with h5py.File(voila_file, 'a') as h:
            index_in_h = 'index' in h

            if remove_index and index_in_h:
                voila_log().info('Removing index from HDF5')
                del h['index']

            return index_in_h

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
        with h5py.File(voila_file, 'a') as h:
            h.create_dataset('index', voila_index.shape, data=voila_index)

    def _create_dtype(self, voila_index):
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

    def _heterogen(self):
        """
        Create index for heterogen analysis type.  This is an odd case as there can be more then one voila file.  The
        current solution is to create a separate voila index hdf5 file in the directory of the first voila file.  Order
        of voila files determined by sorting the filename/location of the each voila file in the config.

        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        voila_dir = Path(config.voila_file).parents[0]
        index_file = voila_dir / 'index.hdf5'

        if index_file.exists() and config.force_index:
            os.remove(index_file)

        if config.force_index or not index_file.exists():
            log.info('Creating index: ' + str(index_file))
            voila_index = []

            with ViewSpliceGraph() as sg, ViewHeterogens() as m:
                for lsv_id in m.lsv_ids():
                    het = m.lsv(lsv_id)

                    gene_id = het.gene_id
                    gene = sg.gene(gene_id)
                    gene_name = gene['name']

                    row = (lsv_id, gene_id, gene_name)

                    lsv_f = [getattr(het, f) for f in lsv_filters]

                    # For some reason, numpy needs these in tuples.
                    row = tuple(chain(row, lsv_f))

                    voila_index.append(row)

            dtype = self._create_dtype(voila_index)
            self._write_index(index_file, voila_index, dtype)

    def _deltapsi(self):
        """
        Generates index for delta psi analysis type.

        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = config.force_index

        for voila_file in config.voila_files:

            if not self._index_in_voila(voila_file, remove_index) or force_index:

                log.info('Creating index: ' + voila_file)
                voila_index = []

                with ViewSpliceGraph() as sg, ViewDeltaPsi() as m:
                    for lsv_id in m.lsv_ids():
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
                        voila_index.append(row)

                dtype = self._create_dtype(voila_index)
                self._write_index(voila_file, voila_index, dtype)

    def _psi(self):
        """
        Create index to PSI analysis type.
        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = config.force_index

        for voila_file in config.voila_files:

            if not self._index_in_voila(voila_file, remove_index) or force_index:

                log.info('Creating index: ' + voila_file)
                voila_index = []

                with ViewSpliceGraph() as sg, ViewPsi() as m:
                    for lsv_id in m.lsv_ids():
                        psi = m.lsv(lsv_id)

                        gene_id = psi.gene_id
                        gene = sg.gene(gene_id)
                        gene_name = gene['name']
                        lsv_f = [getattr(psi, f) for f in lsv_filters]
                        row = (lsv_id, gene_id, gene_name)

                        # for some reason, numpy needs these to tuples.
                        row = tuple(chain(row, lsv_f))

                        voila_index.append(row)

                dtype = self._create_dtype(voila_index)
                self._write_index(voila_file, voila_index, dtype)

    @classmethod
    def psi(cls, gene_id=None):
        """
        Get PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        index_file = ViewConfig().voila_file
        yield from cls.row_data(index_file, gene_id, psi_keys)

    @classmethod
    def delta_psi(cls, gene_id=None):
        """
        Get Delta PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        index_file = ViewConfig().voila_file
        yield from cls.row_data(index_file, gene_id, dpsi_keys)

    @classmethod
    def heterogen(cls, gene_id=None):
        """
        Get Heterogen index data in a dictionary for each row.
        :return:
        """

        voila_file = ViewConfig().voila_file
        voila_file = Path(voila_file)
        voila_dir = voila_file.parents[0]
        index_file = voila_dir / 'index.hdf5'
        yield from cls.row_data(index_file, gene_id, het_keys)

    @staticmethod
    def row_data(index_file, gene_id, keys):
        """
        For each row in index, zip list of keys with values in the row.
        :param gene_id:
        :param keys:
        :return:
        """

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
