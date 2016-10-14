import numpy as np
import scipy.sparse

import majiq.src.config as majiq_config
from majiq.src.constants import *


class Junction:
    __eq__ = lambda self, other: self.start == other.start and self.end == other.end
    __ne__ = lambda self, other: self.start != other.start or self.end != other.end
    __lt__ = lambda self, other: self.start < other.start or (self.start == other.start and self.end < other.end)
    __le__ = lambda self, other: self.start <= other.start or (self.start == other.start and self.end <= other.end)
    __gt__ = lambda self, other: self.start > other.start or (self.start == other.start and self.end > other.end)
    __ge__ = lambda self, other: self.start >= other.start or (self.start == other.start and self.end >= other.end)

    def __init__(self, start, end, donor, acceptor, gene_id, annotated=False, retrieve=False, num_exp=1):
        """ The start and end in junctions are the last exon in """

        self.gene_name = gene_id
        self.id = "%s:%s-%s" % (self.gene_name, start, end)
        self.start = start
        self.end = end
        if donor is None:
            self.donor_id = -1
        else:
            self.donor_id = donor.get_id()
        if acceptor is None:
            self.acceptor_id = -1
        else:
            self.acceptor_id = acceptor.get_id()
        self.annotated = annotated

        if retrieve:
            # self.coverage = scipy.sparse.lil_matrix((num_exp, (majiq_config.readLen - 16) + 1), dtype=np.float)
            self.coverage = scipy.sparse.lil_matrix((num_exp, (majiq_config.readLen - 16) + 1), dtype=np.float)
            self.gc_content = np.zeros((1, (majiq_config.readLen - 16) + 1), dtype=np.float)

        self.transcript_id_list = []

    def __hash__(self):
        return hash(self.start) ^ hash(self.end) ^ hash(self.gene_name)

    def to_db_hdf5(self, hdf5grps):
        h_jnc = hdf5grps.create_group(self.id)
        h_jnc.attrs['id'] = self.id
        if self.start is None:
            h_jnc.attrs['start'] = -1
        else:
            h_jnc.attrs['start'] = self.start

        if self.end is None:
            h_jnc.attrs['end'] = -1
        else:
            h_jnc.attrs['end'] = self.end

        h_jnc.attrs['donor_id'] = self.donor_id
        h_jnc.attrs['acceptor_id'] = self.acceptor_id

    def to_rna_hdf5(self, hdf5grps, dataset, data_index):
        h_jnc = hdf5grps.create_group(self.id)
        h_jnc.attrs['id'] = self.id
        h_jnc.attrs['start'] = self.start
        h_jnc.attrs['end'] = self.end
        h_jnc.attrs['annotated'] = self.annotated

        h_jnc.attrs['chrom'] = self.get_gene().get_chromosome()
        h_jnc.attrs['strand'] = self.get_gene().get_strand()
        h_jnc.attrs['gene_name'] = self.gene_name

        try:
            if data_index == dataset.shape[0]:
                shp_new = dataset.shape[0] + majiq_config.nrandom_junctions
                dataset.resize((shp_new, dataset.shape[1]))

            dataset[data_index, :] = self.coverage[0, :]
            h_jnc.attrs['coverage_index'] = data_index

        except:
            print "HDF5 ERROR", self.id, self.coverage.shape
            raise

    # GETTERs
    def get_id(self):
        return self.id

    def get_coverage(self):
        return self.coverage

    def get_ss_5p(self):
        return self.start

    def get_ss_3p(self):
        return self.end

    def get_gene(self):
        return majiq_config.gene_tlb[self.gene_name]

    def get_coordinates(self):
        return self.start, self.end

    def get_donor(self):
        if self.donor_id == -1:
            ex = None
        elif self.donor_id is None:
            ex = self.get_gene().get_exon_in_coord(self.start)
            if ex is not None:
                self.donor_id = ex.get_id()
        else:
            ex = self.get_gene().get_exon_by_id(self.donor_id)
        return ex

    def is_virtual(self):
        if self.get_acceptor() is None:
            acc_intron = False
        else:
            acc_intron = self.get_acceptor().is_intron()

        if self.get_donor() is None:
            don_intron = False
        else:
            don_intron = self.get_donor().is_intron()

        if acc_intron or don_intron:
            return True
        else:
            return False

    def get_acceptor(self):
        if self.acceptor_id == -1:
            ex = None
        elif self.acceptor_id is None:
            ex = self.get_gene().get_exon_in_coord(self.end)
            if ex is not None:
                self.acceptor_id = ex.get_id()
        else:
            ex = self.get_gene().get_exon_by_id(self.acceptor_id)
        return ex

    def get_gc_content(self):
        return self.gc_content

    def get_read_num(self, idx=0):
        if idx == -1:
            res = self.coverage.sum()
        else:
            res = self.coverage[idx, :].sum()
        return res

    def get_transcript_list(self):
        return self.transcript_id_list

    def is_annotated(self):
        return self.annotated

    def is_reliable(self):
        cov = self.get_coverage()
        res = False
        for tissue, list_idx in majiq_config.tissue_repl.items():
            mu = np.mean(cov[list_idx].sum(axis=1))
            if mu > majiq_config.min_denovo:
                res = True
                break
        return res

    def is_reliable_in_tissue(self):
        return self.get_coverage().sum() > majiq_config.min_denovo

    # MODIFIERs
    def add_gc_content_positions(self, gc):
        self.gc_content = gc

    def add_donor(self, donor):
        if donor is None:
            self.donor_id = -1
        else:
            self.donor_id = donor.get_id()

    def add_acceptor(self, acceptor):
        if acceptor is None:
            self.acceptor_id = -1
        else:
            self.acceptor_id = acceptor.get_id()

    def add_transcript(self, trnscrpt):
        self.transcript_id_list.append(trnscrpt.get_id())

    def update_junction_read(self, read_n, start, gc, unique):
        left_ind = majiq_config.readLen - (self.start - start) - MIN_BP_OVERLAP + 1
        try:
            if unique:
                self.coverage[0, left_ind] += read_n
            else:
                self.coverage[0, left_ind] = -1
            self.gc_content[0, left_ind] = gc
        except:
            print self.gene_name, start, left_ind, self.gc_content.shape
            raise

    def reset_coverage(self):
        self.coverage[:] = 0

    def set_coverage(self, exp_idx, cov):
        self.coverage[exp_idx] = cov