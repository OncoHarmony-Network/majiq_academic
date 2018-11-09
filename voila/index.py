from pathlib import Path

import h5py
import numpy as np

from voila import constants
from voila.api.view_matrix import ViewHeterogens, ViewDeltaPsi, ViewPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.config import Config
from voila.exceptions import UnknownAnalysisType
from voila.utils.voila_log import voila_log


class Index:
    def __init__(self):
        analysis_type = Config().analysis_type

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
        config = Config()
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
        config = Config()
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
                        gene_name = gene.name

                        excl_incl = dpsi.excl_incl
                        excl_incl = max(abs(a - b) for a, b in excl_incl)

                        row = (lsv_id, gene_id, gene_name, excl_incl)
                        self._check_strings(*row)
                        voila_index.append(row)

                voila_index = np.array(voila_index, dtype=np.dtype('S250,S250,S250,f4'))
                with h5py.File(voila_file, 'a') as h:
                    h.create_dataset('index', voila_index.shape, data=voila_index)

    def _psi(self):
        config = Config()
        log = voila_log()
        force_index = remove_index = config.force_index

        for voila_file in config.voila_files:

            if not self._index_in_voila(voila_file, remove_index) or force_index:

                log.info('Creating index: ' + voila_file)
                voila_index = []

                with ViewSpliceGraph() as sg, ViewPsi() as m:
                    for lsv_id in m.lsv_ids():
                        dpsi = m.lsv(lsv_id)

                        gene_id = dpsi.gene_id
                        gene = sg.gene(gene_id)
                        gene_name = gene.name

                        row = (lsv_id, gene_id, gene_name)

                        self._check_strings(*row)

                        voila_index.append(row)

                voila_index = np.array(voila_index, dtype=np.dtype('S250,S250,S250'))
                with h5py.File(voila_file, 'a') as h:
                    h.create_dataset('index', voila_index.shape, data=voila_index)
