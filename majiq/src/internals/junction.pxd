from libcpp.string cimport string

cdef extern from "junction.hpp":
    cdef cppclass Junction:
        Junction() except +
        Junction(string chrom1, char strand1, unsigned int start1, unsigned int end1) except +
        unsigned int nreads;
        int index;
        int lsv_index;
        bint annot;
        string chrom;
        char strand;
        unsigned int start;
        unsigned int end;
