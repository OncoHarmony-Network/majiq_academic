from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.list cimport list

cdef extern from "grimoire.hpp" namespace "grimoire":
    cdef cppclass Junction:
        Junction() except +
        Junction(string gene_id1, int start1, int end1, int nexp, int eff_len) except +
        # Junction(string gene_id1, int start1, int end1) except +
        unsigned int nreads;
        string gene_id;
        unsigned int n_exps;
        int start;
        int end;

    cdef cppclass Gene:
        Gene() except +
        Gene(string id1, string name1, string chromosome1,
             char strand1, unsigned int start1, unsigned int end1) except +
        string id;
        string name;
        string chromosome;
        char strand;
        unsigned int start;
        unsigned int end;

    cdef cppclass Exon:
        Exon() except +
        Exon(int start1, int end1) except +
        Exon(int start1, int end1, bint annot1) except +
        Exon(int start1, int end1, bint annot1, bint intron1) except +

        int start ;
        int end ;
        bint annot ;
        unsigned int db_start ;
        unsigned int db_end ;
        set[Junction *] ib ;
        set[Junction *] ob ;
        bint intron ;


    cdef cppclass LSV:
        string id ;
        set[Junction *] junctions ;
        string type ;
        LSV() except +
        LSV(string gene_id1, char strand, Exon* ex, bint ss) except +


    int detect_lsvs(list[LSV*] out_lsvlist, map[string, Exon*] exon_map, Gene* gObj, unsigned int nexp, unsigned int eff_len,
                    int minpos, int minreads) nogil ;
    void detect_exons(map[string, Junction*] junc_map, map[string, Exon*] exon_map) nogil ;
    void free_gene(Gene * gObj, map[string, Junction*] junc_map, map[string, Exon*] exon_map) nogil ;
    float * boostrap_samples(LSV* lsvObj, int msamples, int ksamples, int exp_idx, int eff_len) nogil ;