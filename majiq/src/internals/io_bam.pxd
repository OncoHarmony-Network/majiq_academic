from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.set cimport set
from majiq.src.internals.junction cimport Junction

cdef extern from "io_bam.hpp" namespace "io_bam":
    cdef cppclass IOBam:
        IOBam() except +
        IOBam(char *, int) except +
        IOBam(char*, int, char*) except +
        int find_junctions()
        map[string, Junction] get_dict()
        void set_filters(set[string], set[string])
        map[string, Junction] junc_dict_


