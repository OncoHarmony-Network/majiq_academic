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


    cdef cppclass Gene:
        Gene() nogil except +
        Gene(string id1, string name1, string chromosome1,
             char strand1, unsigned int start1, unsigned int end1) nogil except +
        # ~Gene() except +

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
                                  unsigned int denovo_thresh, unsigned int min_experiments) nogil ;
        void    fill_junc_tlb(map[string, vector[string]]& tlb) nogil ;
        void    connect_introns() nogil ;

        map[string, Junction*] junc_map_ ;
        map[string, Exon*] exon_map_ ;
        vector[Intron*] intron_vec_ ;

    cdef cppclass Exon:
        Exon() nogil except +
        Exon(int start1, int end1) nogil except +
        Exon(int start1, int end1, bint annot1) nogil except +
        int     get_start() nogil ;
        int     get_end()   nogil ;
        bint annot_ ;
        int db_start_ ;
        int db_end_ ;
        set[Junction *] ib ;
        set[Junction *] ob ;
        Intron * ob_irptr ;

        bint has_out_intron() nogil ;



    cdef cppclass LSV:
        set[Junction *] junctions ;
        LSV() nogil except +
        LSV(string gene_id1, char strand, Exon* ex, bint ss) nogil except +
        bint    gather_lsv_info(float* source, float* target, clist[Jinfo*]& info, map[string, Jinfo]& tlb,
                             unsigned int msample) nogil ;
        string  get_id() nogil ;
        Gene*   get_gene() nogil ;
        string  get_type() nogil ;
        # int    get_num_junctions() nogil ;
        int     get_num_variations() nogil ;
        vector[Junction *] get_junctions() nogil ;
        Intron * get_intron() nogil ;

    cdef struct Jinfo:
        unsigned int index ;
        unsigned int start ;
        unsigned int end ;
        int sreads ;
        int npos ;

    int detect_lsvs(vector[LSV*] out_lsvlist, Gene* gObj) nogil ;
    float * boostrap_samples(LSV* lsvObj, int msamples, int ksamples, int exp_idx, int eff_len) nogil ;
    void sortGeneList(vector[Gene*] glist) nogil ;
    vector[Intron *]  find_intron_retention(vector[Gene*]& gene_list, int start, int end) nogil ;



