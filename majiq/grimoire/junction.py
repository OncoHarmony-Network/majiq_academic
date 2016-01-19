import numpy as np
import scipy.sparse

import majiq.src.config as majiq_config


class Junction:
    __eq__ = lambda self, other: self.start == other.start and self.end == other.end
    __ne__ = lambda self, other: self.start != other.start or self.end != other.end
    __lt__ = lambda self, other: self.start < other.start or (self.start == other.start and self.end < other.end)
    __le__ = lambda self, other: self.start <= other.start or (self.start == other.start and self.end <= other.end)
    __gt__ = lambda self, other: self.start > other.start or (self.start == other.start and self.end > other.end)
    __ge__ = lambda self, other: self.start >= other.start or (self.start == other.start and self.end >= other.end)

    def __init__(self, start, end, donor, acceptor, gene, readN=0, annotated=False):
        ''' The start and end in junctions are the last exon in '''

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
        self.gene_name = gene.get_id()
        self.annotated = annotated
        self.coverage = scipy.sparse.lil_matrix((majiq_config.num_experiments,  (majiq_config.readLen - 16) + 1),
                                                dtype=np.float)
        self.gc_content = scipy.sparse.lil_matrix((1, (majiq_config.readLen - 16) + 1), dtype=np.float)
        self.id = "%s:%s-%s" % (self.gene_name, start, end)
        self.transcript_id_list = []

    def __hash__(self):
        return hash(self.start) ^ hash(self.end) ^ hash(self.gene_name)

    # def __del__(self):
    # junc_list = self.get_gene().get_junction_list()

    # GETTERs
    def get_id(self):
        return self.id

    def get_coverage(self, experiment=None):
        if experiment is None:
            return self.coverage
        else:
            return self.coverage[experiment, :]

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
            if not ex is None:
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
            if not ex is None:
                self.acceptor_id = ex.get_id()
        else:
            ex = self.get_gene().get_exon_by_id(self.acceptor_id)
        return ex

    def get_gc_content(self):
        return self.gc_content

    def get_read_num(self, idx):
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
        cov = self.get_coverage().toarray()
        res = False
        for tissue, list_idx in majiq_config.tissue_repl.items():
            mu = np.mean(cov[list_idx].sum(axis=1))
            if mu > majiq_config.min_denovo:
                res = True
                break
        return res

    # MODIFIERs

    def add_gc_content_positions(self, pos, gc):
        self.gc_content = gc

    #        self.gc_factor[exp_idx,:] = gc_factor

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

    def update_junction_read(self, exp_idx, read_n, start, gc, unique):
        left_ind = majiq_config.readLen - (self.start - start) - 8 + 1
        if unique:
            self.coverage[exp_idx, left_ind] += read_n
        else:
            self.coverage[exp_idx, left_ind] = -1
        self.gc_content[0, left_ind] = gc

    def reset_coverage(self, exp_idx):
        self.coverage[exp_idx, :] = 0


    def hdf5_junction(self, hdf5grps, exp_idx):
        h_jnc = hdf5grps.create_group(self.get_id())
        self.annot = self.is_annotated()
        if self is None:
            #self.gc_index = scipy.sparse.lil_matrix((1, (max(config.readLen) - 16) + 1), dtype=np.int)
            self.name = None
            self.id = "None"

            h_jnc['chrom'] = None
            h_jnc['coord1'] = [0, 0]
            h_jnc['coord2'] = [0, 0]
            h_jnc['strand'] = None
            h_jnc['coverage'] = scipy.sparse.lil_matrix((1, (majiq_config.readLen - 16) + 1), dtype=np.float)
        else:
            h_jnc['name'] = "%s:%s-%s" % (self.get_gene().get_id(), self.get_ss_5p(), self.get_ss_3p())
            h_jnc['id'] = "%s:%s-%s" % (self.get_gene().get_chromosome(), self.get_ss_5p(), self.get_ss_3p())

            h_jnc['chrom'] = self.get_gene().get_chromosome()
            h_jnc['strand'] = self.get_gene().get_strand()
            if self.get_donor() is None:
                h_jnc['coord1'] = [0, 0]
            else:
                h_jnc['coord1'] = self.get_donor().get_coordinates()

            if self.get_acceptor() is None:
                h_jnc['coord2'] = [0, 0]
            else:
                h_jnc['coord2'] = self.get_acceptor().get_coordinates()

            h_jnc.create_dataset('coverage', data=self.coverage[exp_idx, :].toarray())
            if majiq_config.gcnorm:
                gc_factor = scipy.sparse.lil_matrix((1, (majiq_config.readLen - 16) + 1), dtype=np.float)
            else:
                gc_factor = None
            for jj in range(majiq_config.readLen - 16 + 1):
                dummy = self.gc_content[0, jj]
                gc_factor[0, jj] *= dummy
            h_jnc.create_dataset('gc_factor', data=gc_factor.toarray())

        return h_jnc


def set_gc_factor(hdf5grp, exp_idx):
    if majiq_config.gcnorm:
        gc_factor = hdf5grp['gc_factor']
        nnz = gc_factor.nonzero()
        for idx in xrange(nnz[0].shape[0]):
            i = nnz[0][idx]
            j = nnz[1][idx]
            dummy = gc_factor[i, j]
            hdf5grp['coverage'][i, j] *= majiq_config.gc_factor[exp_idx](dummy)
    del hdf5grp['gc_factor']


class MajiqJunction:
    def __init__(self, jnc, exp_idx):
        self.exons = {}
        self.annot = jnc.is_annotated()
        if jnc is None:
            #self.gc_index = scipy.sparse.lil_matrix((1, (max(config.readLen) - 16) + 1), dtype=np.int)
            self.name = None
            self.id = "None"

            self.exons['chrom'] = None
            self.exons['coord1'] = [0, 0]
            self.exons['coord2'] = [0, 0]
            self.exons['strand'] = None
            self.coverage = scipy.sparse.lil_matrix((1, (majiq_config.readLen - 16) + 1), dtype=np.float)
        else:
            self.name = "%s:%s-%s" % (jnc.get_gene().get_id(), jnc.get_ss_5p(), jnc.get_ss_3p())
            self.id = "%s:%s-%s" % (jnc.get_gene().get_chromosome(), jnc.get_ss_5p(), jnc.get_ss_3p())

            self.exons['chrom'] = jnc.get_gene().get_chromosome()
            self.exons['strand'] = jnc.get_gene().get_strand()
            if jnc.get_donor() is None:
                self.exons['coord1'] = [0, 0]
            else:
                self.exons['coord1'] = jnc.get_donor().get_coordinates()

            if jnc.get_acceptor() is None:
                self.exons['coord2'] = [0, 0]
            else:
                self.exons['coord2'] = jnc.get_acceptor().get_coordinates()

            self.coverage = jnc.coverage[exp_idx, :]
            if majiq_config.gcnorm:
                self.gc_factor = scipy.sparse.lil_matrix((1, (majiq_config.readLen - 16) + 1), dtype=np.float)
            else:
                self.gc_factor = None
            for jj in range(majiq_config.readLen - 16 + 1):
                dummy = jnc.gc_content[0, jj]
                self.gc_factor[0, jj] *= dummy

    def set_gc_factor(self, exp_idx):
        if majiq_config.gcnorm:
            nnz = self.gc_factor.nonzero()
            for idx in xrange(nnz[0].shape[0]):
                i = nnz[0][idx]
                j = nnz[1][idx]
                dummy = self.gc_factor[i, j]
                self.coverage[i, j] *= majiq_config.gc_factor[exp_idx](dummy)
        del self.gc_factor


