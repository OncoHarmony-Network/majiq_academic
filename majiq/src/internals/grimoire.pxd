from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.list cimport list as clist
from majiq.src.internals.psi cimport psi_distr_t

ctypedef vector[overGene*] overGene_vect_t
ctypedef vector[Gene*] Gene_vect_t
cdef extern from "grimoire.hpp" namespace "grimoire":
    cdef cppclass overGene:
        pass

    cdef cppclass Junction:
        Junction() nogil except +
        Junction(int start1, int end1) nogil except +
        Junction(int start1, int end1, bint annot1) nogil except +
        void         update_flags(unsigned int num_reads, unsigned int num_pos, unsigned int denovo_thresh,
                                  unsigned int min_experiments, bint denovo) nogil
        void         clear_nreads(bint reset_grp) nogil
        int          get_start() nogil ;
        int          get_end() nogil ;
        string       get_key() nogil ;
        string       get_key(Gene * gObj) nogil ;
        string       get_key(Gene * gObj, int strandness) nogil ;
        bint         get_intronic() nogil
        bint         get_annot() nogil
        bint         get_bld_fltr() nogil
        bint         get_denovo_bl() nogil

        unsigned int nreads;


    cdef cppclass Intron:

            Intron () nogil except +
            Intron (int start1, int end1, bint annot1, Gene* gObj1) nogil except +
            int     get_start()  nogil ;
            int     get_end()    nogil ;
            bint    get_annot()  nogil ;
            string  get_key() nogil ;
            string  get_key(Gene * gObj) nogil ;
            Gene *  get_gene() nogil ;
            bint    get_ir_flag() nogil ;


    cdef cppclass Gene:
        Gene() nogil except +
        Gene(string id1, string name1, string chromosome1,
             char strand1, unsigned int start1, unsigned int end1) nogil except +

        string  get_chromosome() nogil ;
        int     get_start()      nogil ;
        int     get_strand()     nogil ;
        int     get_end()        nogil ;
        string  get_id()         nogil ;
        string  get_name()       nogil ;
        void    print_gene()     nogil ;
        void    print_exons()    nogil ;
        void    detect_exons()   nogil ;
        void    detect_introns() nogil ;

        void    create_annot_intron(int start_ir, int end_ir) nogil ;
        void    add_elements(map[string, Junction*] junc_map, map[string, Exon*] exon_map) nogil ;
        void    update_junc_flags(int efflen, bint is_last_exp, unsigned int minreads, unsigned int minpos,
                                  unsigned int denovo_thresh, unsigned int min_experiments, bint denovo) nogil ;
        void    fill_junc_tlb(map[string, vector[string]]& tlb) nogil ;
        void    connect_introns() nogil ;
        int     detect_lsvs(vector[LSV*] out_lsvlist) nogil ;

        map[string, Junction*] junc_map_ ;
        map[string, Exon*] exon_map_ ;
        vector[Intron*] intron_vec_ ;


    cdef cppclass Exon:
        Exon() nogil except +
        Exon(int start1, int end1) nogil except +
        Exon(int start1, int end1, bint annot1) nogil except +
        int     get_start() nogil ;
        int     get_end()   nogil ;
        bint    has_out_intron() nogil ;
        bint annot_ ;
        int db_start_ ;
        int db_end_ ;
        set[Junction *] ib ;
        set[Junction *] ob ;
        Intron * ob_irptr ;


    cdef cppclass LSV:
        set[Junction *] junctions ;
        LSV() nogil except +
        LSV(string gene_id1, char strand, Exon* ex, bint ss) nogil except +
        bint                gather_lsv_info(float* source, float* target, clist[Jinfo*]& info, map[string, Jinfo]& tlb,
                                            unsigned int msample) nogil ;
        string              get_id() nogil ;
        Gene*               get_gene() nogil ;
        string              get_type() nogil ;
        int                 get_num_variations() nogil ;
        vector[Junction *]  get_junctions() nogil ;
        Intron *            get_intron() nogil ;


    cdef cppclass Jinfo:
        unsigned int index ;
        int sreads ;
        int npos ;
        Jinfo() nogil ;
        Jinfo(unsigned int index1, int sreads1, int npos1) nogil ;



    # cdef cppclass qLSV:
    #     vector[vector[psi_distr_t]] coverages ;
    #
    #     qLSV() nogil ;
    #     qLSV(vector[string] & lsv_list_) nogil ;
    #         #~qLSV()                        { omp_destroy_lock(&semaphore_lock) ; }
    #     vector[string] & get_lsvlist() nogil ;
    #     void free_lock() nogil ;
    #     void test_lock() nogil ;
    #     void add_counter() nogil ;



    void sortGeneList(vector[Gene*] glist) nogil ;
    # vector[Intron *]  find_intron_retention(vector[Gene*]& gene_list, char strand, int start, int end) nogil ;
    # vector[Intron *]  find_intron_retention(vector[Gene*]& gene_list, string geneid, int start, int end) nogil ;
    vector[Intron *]  find_intron_retention(Gene * gObj, int start, int end) nogil ;
    # vector[Gene*] find_gene_from_junc(const map[string, overGene_vect_t]& glist, string chrom, int start, int end, bint ir) nogil ;
    void find_gene_from_junc(map[string, overGene_vect_t] & glist, string chrom, char strand, int start, int end,
                             vector[Gene*]& oGeneList, bint ir) nogil ;
    void fill_junc_tlb(vector[LSV*]& lsv_list, map[string, int]& tlb) nogil ;
    bint isNullJinfo(Jinfo* x) nogil ;
    void free_JinfoVec( vector[Jinfo*]& jvec) nogil ;
    string key_format(string gid, int coord1, int coord2, bint ir) nogil ;

