import json
from pathlib import Path

import h5py
import numpy as np

from voila import constants
from voila.api.view_matrix import ViewHeterogens, ViewDeltaPsi, ViewPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.config import ViewConfig
from voila.exceptions import UnknownAnalysisType, IndexNotFound
from voila.utils.voila_log import voila_log
from voila.vlsv import matrix_area

lsv_filters = ['a5ss', 'a3ss', 'exon_skipping', 'target', 'source', 'binary', 'complex']
psi_keys = ['lsv_id', 'gene_id', 'gene_name'] + lsv_filters
dpsi_keys = ['lsv_id', 'gene_id', 'gene_name', 'excl_incl', 'dpsi_threshold', 'confidence_threshold'] + lsv_filters


class Index:
    def __init__(self):
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
    def _check_strings(*strings):
        assert any(len(str(s)) <= 250 for s in strings)

    @staticmethod
    def _index_in_voila(voila_file, remove_index=False):
        with h5py.File(voila_file, 'a') as h:
            index_in_h = 'index' in h

            if remove_index and index_in_h:
                del h['index']
                voila_log().info('Removing index from HDF5')

            return index_in_h

    def _heterogen(self):
        config = ViewConfig()
        log = voila_log()

        voila_dir = Path(config.voila_file).parents[0]
        index_file = voila_dir / 'index.hdf5'

        if config.force_index or not index_file.exists():
            log.info('Creating index: ' + str(index_file))
            voila_index = []

            with ViewSpliceGraph() as sg, ViewHeterogens() as m:
                for lsv_id in m.lsv_ids():
                    het = m.lsv(lsv_id)

                    gene_id = het.gene_id
                    gene = sg.gene(gene_id)
                    gene_name = gene.name

                    row = (lsv_id, gene_id, gene_name)
                    self._check_strings(*row)

                    voila_index.append(row)

            voila_index = np.array(voila_index, dtype=np.dtype('S250,S250,S250'))
            with h5py.File(index_file, 'w') as h:
                h.create_dataset('index', voila_index.shape, data=voila_index, compression='gzip', compression_opts=9)

            log.info('Index created')

    def _deltapsi(self):
        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = config.force_index

        for voila_file in config.voila_files:

            if not self._index_in_voila(voila_file, remove_index) or force_index:

                log.info('Creating index: ' + voila_file)
                voila_index = []
                lsv_id_len = 0
                gene_id_len = 0
                gene_name_len = 0
                dpsi_thresh_len = 0
                confidence_thresh_len = 0

                with ViewSpliceGraph() as sg, ViewDeltaPsi() as m:
                    for lsv_id in m.lsv_ids():
                        dpsi = m.lsv(lsv_id)

                        gene_id = dpsi.gene_id
                        gene = sg.gene(gene_id)
                        gene_name = gene.name

                        dpsi_thresh = dpsi.means
                        dpsi_thresh = np.abs(dpsi_thresh)
                        dpsi_thresh = dpsi_thresh.tolist()
                        dpsi_thresh = json.dumps(dpsi_thresh)

                        bins = dpsi.bins
                        confidence_thresh = list(max(matrix_area(b, x) for b in bins) for x in np.linspace(0, 1, 10))
                        confidence_thresh = json.dumps(confidence_thresh)

                        lsv_id_len = max(lsv_id_len, len(lsv_id))
                        gene_id_len = max(gene_id_len, len(gene_id))
                        gene_name_len = max(gene_name_len, len(gene_name))
                        dpsi_thresh_len = max(dpsi_thresh_len, len(dpsi_thresh))
                        confidence_thresh_len = max(confidence_thresh_len, len(confidence_thresh))

                        excl_incl = dpsi.excl_incl
                        excl_incl = max(abs(a - b) for a, b in excl_incl)

                        lsv_f = [getattr(dpsi, f) for f in lsv_filters]

                        row = tuple([lsv_id, gene_id, gene_name, excl_incl, dpsi_thresh, confidence_thresh] + lsv_f)
                        voila_index.append(row)

                dtype = 'S{},S{},S{},f4,S{},S{},?,?,?,?,?,?,?'.format(lsv_id_len, gene_id_len, gene_name_len,
                                                                      dpsi_thresh_len, confidence_thresh_len)
                voila_index = np.array(voila_index, dtype=np.dtype(dtype))
                with h5py.File(voila_file, 'a') as h:
                    h.create_dataset('index', voila_index.shape, data=voila_index)

    def _psi(self):
        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = config.force_index

        for voila_file in config.voila_files:

            if not self._index_in_voila(voila_file, remove_index) or force_index:

                log.info('Creating index: ' + voila_file)
                voila_index = []
                lsv_id_len = 0
                gene_id_len = 0
                gene_name_len = 0

                with ViewSpliceGraph() as sg, ViewPsi() as m:
                    for lsv_id in m.lsv_ids():
                        psi = m.lsv(lsv_id)

                        gene_id = psi.gene_id
                        gene = sg.gene(gene_id)
                        gene_name = gene.name
                        lsv_f = [getattr(psi, f) for f in lsv_filters]

                        lsv_id_len = max(lsv_id_len, len(lsv_id))
                        gene_id_len = max(gene_id_len, len(gene_id))
                        gene_name_len = max(gene_name_len, len(gene_name))

                        # for some reason, numpy needs these to tuples.
                        row = tuple([lsv_id, gene_id, gene_name] + lsv_f)

                        voila_index.append(row)

                dtype = 'S{},S{},S{},?,?,?,?,?,?,?'.format(lsv_id_len, gene_id_len, gene_name_len)
                voila_index = np.array(voila_index, dtype=np.dtype(dtype))

                with h5py.File(voila_file, 'a') as h:
                    h.create_dataset('index', voila_index.shape, data=voila_index)

    @classmethod
    def psi(cls, gene_id=None):
        yield from cls.row_data(gene_id, psi_keys)

    @classmethod
    def delta_psi(cls, gene_id=None):
        yield from cls.row_data(gene_id, dpsi_keys)

    @staticmethod
    def row_data(gene_id, keys):
        try:
            gene_id = gene_id.encode('utf-8')
        except AttributeError:
            pass

        config = ViewConfig()

        with h5py.File(config.voila_file, 'r') as h:
            try:
                for row in h['index'].value:
                    if gene_id is None or gene_id == row[1]:
                        yield dict(zip(keys, row))
            except KeyError:
                raise IndexNotFound()
