from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector
from majiq.src.internals.grimoire cimport Junction, Gene

cdef extern from "io_bam.hpp" namespace "io_bam":
    cdef cppclass IOBam:
        IOBam() nogil except +
        IOBam(string, int, unsigned int, unsigned int, unsigned int) nogil except +
        IOBam(string, int, unsigned int, string) nogil except +
        int find_junctions(int min_experiments) nogil
        # int find_junctions_from_region(Gene * gobj) nogil
        int find_junctions_from_region(vector[Gene*] gobj) nogil
        map[string, Junction] get_dict() nogil
        void set_filters(map[string, set[string]]*) nogil
        int ParseJunctionsFromFile(string filename, int nthreads) nogil

