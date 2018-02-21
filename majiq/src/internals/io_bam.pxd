from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.set cimport set
from majiq.src.internals.junction cimport Junction

cdef extern from "io_bam.hpp" namespace "io_bam":
    cdef cppclass IOBam:
        IOBam() nogil except +
        IOBam(string, int, map[string, Junction]*) nogil except +
        IOBam(string, int, map[string, Junction]*, char*) nogil except +
        int find_junctions(int min_experiments) nogil
        map[string, Junction] get_dict() nogil
        void set_filters(map[string, set[string]]*) nogil
        map[string, Junction] junc_dict_

