from libcpp.string cimport string

cdef extern from "junction.hpp":
    cdef cppclass Junction:
        Junction() except +
        Junction(string gene_id1, unsigned int start1, unsigned int end1) except +
        unsigned int nreads;
        string gene_id;
        unsigned int n_exps;
        unsigned int start;
        unsigned int end;
