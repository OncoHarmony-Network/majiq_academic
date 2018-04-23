from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.list cimport list as clist

cdef extern from "grimoire.hpp" namespace "grimoire":
    cdef cppclass Junction:
        Junction() nogil except +
        Junction(int start1, int end1) nogil except +
        Junction(int start1, int end1, bint annot1) nogil except +
        void         update_flags(unsigned int num_reads, unsigned int num_pos, unsigned int denovo_thresh,
                                  unsigned int min_experiments) nogil
        void         clear_nreads(bint reset_grp) nogil
        int          get_start() nogil ;
        int          get_end() nogil ;
        string       get_key() nogil ;
        string       get_key(Gene * gObj) nogil ;
        unsigned int get_sum() nogil
        unsigned int get_npos() nogil
        bint         get_intronic() nogil
        bint         get_annot() nogil
        bint         get_bld_fltr() nogil

        # ~Junction() except +
        # Junction(string gene_id1, int start1, int end1) except +
        unsigned int nreads;
        int start;
        int end;

    cdef cppclass Gene:
        Gene() nogil except +
        Gene(string id1, string name1, string chromosome1,
             char strand1, unsigned int start1, unsigned int end1) nogil except +
        # ~Gene() except +

        void add_elements(map[string, Junction*] junc_map, map[string, Exon*] exon_map) nogil ;
        void detect_exons() nogil ;
        void print_gene() nogil ;
        string get_chromosome() nogil ;
        int get_start() nogil ;
        int get_strand() nogil ;
        int get_end() nogil ;
        string get_id() nogil ;
        string get_name() nogil ;
        void print_exons() nogil ;

        map[string, Junction*] junc_map_ ;
        map[string, Exon*] exon_map_ ;

    cdef cppclass Exon:
        Exon() nogil except +
        Exon(int start1, int end1) nogil except +
        Exon(int start1, int end1, bint annot1) nogil except +
        Exon(int start1, int end1, bint annot1, bint intron1) nogil except +
        int start_ ;
        int end_ ;
        bint annot_ ;
        bint intron_ ;
        int db_start_ ;
        int db_end_ ;
        set[Junction *] ib ;
        set[Junction *] ob ;



    cdef cppclass LSV:
        set[Junction *] junctions ;
        LSV() nogil except +
        LSV(string gene_id1, char strand, Exon* ex, bint ss) nogil except +
        bint gather_lsv_info(float* source, float* target, clist[Jinfo*]& info, map[string, Jinfo]& tlb,
                             unsigned int msample) nogil ;
        string get_id() nogil ;
        Gene*  get_gene() nogil ;
        string get_type() nogil ;
        int    get_num_junctions() nogil ;
        vector[Junction *] get_junctions() nogil ;

    cdef struct Jinfo:
        unsigned int index ;
        unsigned int start ;
        unsigned int end ;
        int sreads ;
        int npos ;

    int detect_lsvs(vector[LSV*] out_lsvlist, Gene* gObj) nogil ;
    float * boostrap_samples(LSV* lsvObj, int msamples, int ksamples, int exp_idx, int eff_len) nogil ;
    void sortGeneList(vector[Gene*] glist) nogil ;



