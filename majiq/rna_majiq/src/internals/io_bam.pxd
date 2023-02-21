from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
cimport numpy as np
from rna_majiq.src.internals.grimoire cimport Junction, Gene, overGene

ctypedef np.float64_t DTYPE_t
ctypedef vector[overGene*] overGene_vect_t

cdef extern from "io_bam.hpp" namespace "io_bam":

    cdef cppclass IOBam:
        IOBam() nogil except +
        IOBam(string, int, unsigned int, unsigned int, map[string, overGene_vect_t], bint simpl1, bint allow_full_intergene) nogil except +
        int ParseJunctionsFromFile(bint ir_func) nogil
        void EstimateEffLenFromFile(int num_reads) nogil
        int bootstrap_samples(int msamples, np.float32_t* boots, np.float32_t pvalue_limit) nogil
        void detect_introns(np.float32_t min_intron_cov, unsigned int min_experiments, np.float32_t min_bins, bint reset,
                            bint denovo) nogil
        int get_njuncs() nogil
        map[string, unsigned int] get_junc_map() nogil
        vector[Junction *] get_junc_vec() nogil
        int * get_junc_vec_summary() nogil ;
        unsigned int get_junc_limit_index() nogil;
        void simplify() nogil ;
        void free_iobam() nogil ;
        void get_intron_raw_cov(np.float32_t* out_cov) nogil ;
        void get_junction_raw_cov(np.float32_t* out_cov) nogil;
        void parseJuncEntry(map[string, overGene_vect_t]& glist, string gid, string chrom, char strand,
                               int start, int end, unsigned int sreads, unsigned int minreads_t, unsigned int numpos,
                               unsigned int minpos_t, unsigned int denovo_t, bint denovo, bint denovo_ir, vector[Gene*]& oGeneList,
                               bint ir, vector[np.float32_t]& ircov, np.float32_t min_intron_cov, np.float32_t min_bins,
                               int minexp, bint reset) nogil ;
        unsigned int get_eff_len() nogil

        vector[shared_ptr[vector[np.float32_t]]] junc_vec ;

    void free_genelist(map[string, overGene_vect_t]& geneList) nogil ;
    void prepare_genelist(map[string, Gene*]& gene_map, map[string, overGene_vect_t]& geneList) nogil ;
    unsigned int eff_len_from_read_length(int read_length) nogil
