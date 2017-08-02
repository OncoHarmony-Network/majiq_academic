import numpy as np
from majiq.src.config import Config
from majiq.src.constants import *


class Junction:
    __eq__ = lambda self, other: self.start == other.start and self.end == other.end
    __ne__ = lambda self, other: self.start != other.start or self.end != other.end
    __lt__ = lambda self, other: self.start < other.start or (self.start == other.start and self.end < other.end)
    __le__ = lambda self, other: self.start <= other.start or (self.start == other.start and self.end <= other.end)
    __gt__ = lambda self, other: self.start > other.start or (self.start == other.start and self.end > other.end)
    __ge__ = lambda self, other: self.start >= other.start or (self.start == other.start and self.end >= other.end)

    def __init__(self, start, end, donor, acceptor, gene_id, annotated=False, retrieve=False, num_exp=1, jindex=-1,
                 intronic=False):
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
            if num_exp == 0:
                num_exp = 1
                majiq_config = Config()
                self.coverage = np.zeros((num_exp, (majiq_config.readLen - 16) + 1), dtype=np.float)
                self.gc_content = np.zeros((1, (majiq_config.readLen - 16) + 1), dtype=np.float)
                self.all_data = True
            else:
                self.idx = jindex
                self.all_data = False
            self.intronic = intronic
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
        h_jnc.attrs['transcript_id_list'] = [xx.encode('utf8') for xx in self.transcript_id_list]

    def to_rna_hdf5(self, hdf5grps, dataset, data_index, gc_dataset=None):
        h_jnc = hdf5grps.create_group(self.id)
        h_jnc.attrs['id'] = self.id
        h_jnc.attrs['start'] = self.start
        h_jnc.attrs['end'] = self.end
        h_jnc.attrs['annotated'] = self.annotated

        h_jnc.attrs['chrom'] = self.get_gene().get_chromosome()
        h_jnc.attrs['strand'] = self.get_gene().get_strand()
        h_jnc.attrs['gene_name'] = self.gene_name

        try:

            dataset.append(self.coverage[0, :])
            if gc_dataset is not None:
                gc_dataset.append(self.gc_content[0, :])

            h_jnc.attrs['coverage_index'] = data_index

        except:
            print("HDF5 ERROR", self.id, self.coverage.shape)
            raise

    # GETTERs
    def get_id(self):
        return self.id

    def get_coverage(self):
        if self.all_data:
            return self.coverage
        else:
            if self.idx == -1:
                majiq_config = Config()
                cov = np.zeros(shape=(majiq_config.num_experiments, 2))
            else:
                cov = self.get_gene().junc_matrix[self.idx]
            return cov

    def get_gene(self):
        majiq_config = Config()
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

    def get_gc_content(self,exp_idx):
        if self.all_data:
            return self.gc_content
        else:
            if self.idx == -1:
                gc = None
            else:
                gc = self.get_gene().gc_content[self.idx, exp_idx]
            return gc

    def get_read_num(self, idx=0):
        if self.all_data:
            cov = self.coverage
        else:
            cov = self.get_gene().junc_matrix[self.idx]

        if idx == -1:
            res = cov.sum(dtype=np.uint32)
        else:
            res = cov[idx].sum(dtype=np.uint32)

        return res

    def get_index(self):
        try:
            return self.idx
        except:
            pass

    def get_gene_name(self):
        return self.gene_name

    def get_coverage_sum(self, idx):

        return self.get_coverage()[idx].sum()

    def get_transcript_list(self):
        return self.transcript_id_list

    def is_reliable(self):
        res = False
        cov = self.get_coverage().sum(axis=1)
        majiq_config = Config()
        for tissue, list_idx in majiq_config.tissue_repl.items():
            mu = np.mean(cov[list_idx])
            if mu > majiq_config.min_denovo:
                res = True
                break
        return res

    def is_reliable_in_tissue(self):
        majiq_config = Config()
        return self.get_coverage().sum() > majiq_config.min_denovo


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
        majiq_config = Config()
        left_ind = majiq_config.readLen - (self.start - start) - MIN_BP_OVERLAP + 1
        try:
            if unique:
                self.coverage[0, left_ind] += read_n
            else:
                self.coverage[0, left_ind] = -1
            self.gc_content[0, left_ind] = gc
        except:
            print(self.gene_name, start, left_ind, self.gc_content.shape)
            raise
