from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
cimport numpy as np
from majiq.src.internals.grimoire cimport Junction, Gene, overGene

ctypedef np.float64_t DTYPE_t
# ctypedef vector[pair[int, vector[Gene*]]] gene_vect_t
# ctypedef vector[Intron] intron_vect_t
ctypedef vector[overGene*] overGene_vect_t

cdef extern from "io_bam.hpp" namespace "io_bam":

    cdef cppclass IOBam:
        IOBam() nogil except +
        IOBam(string, int, unsigned int, unsigned int, map[string, overGene_vect_t], bint simpl1) nogil except +
        # IOBam(string, int, unsigned int) nogil except +
        # int find_junctions(int min_experiments) nogil
        # int find_junctions_from_region(vector[Gene*] gobj) nogil
        int ParseJunctionsFromFile(bint ir_func) nogil
        # int ParseJunctionsFromFile(map[string, ITNode*]) nogil
        int boostrap_samples(int msamples, int ksamples, np.float32_t* boots, float fitfunc_r, float pvalue_limit) nogil
        void detect_introns(float min_intron_cov, unsigned int min_experiments, float min_bins, bint reset) nogil
        int get_njuncs() nogil
        map[string, unsigned int] get_junc_map() nogil
        vector[Junction *] get_junc_vec() nogil
        int * get_junc_vec_summary() nogil ;
        unsigned int get_junc_limit_index() nogil;
        void simplify() nogil ;
        void free_iobam() nogil ;
        void get_intron_raw_cov(np.float32_t* out_cov) nogil ;
        void parseJuncEntry(map[string, overGene_vect_t]& glist, string chrom, char strand, int start, int end,
                            int sreads, vector[Gene*]& oGeneList, bint ir, vector[float]& ircov,
                            float min_intron_cov, float min_bins, int minexp, bint reset) nogil ;

        vector[np.float32_t *] junc_vec ;

    void free_genelist(map[string, overGene_vect_t]& geneList) nogil ;
    void prepare_genelist(map[string, Gene*]& gene_map, map[string, overGene_vect_t]& geneList) nogil ;
