from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector
cimport numpy as np
from majiq.src.internals.grimoire cimport Junction, Gene

ctypedef np.float64_t DTYPE_t

cdef extern from "io_bam.hpp" namespace "io_bam":
    cdef cppclass IOBam:
        IOBam() nogil except +
        IOBam(string, int, unsigned int, unsigned int) nogil except +
        IOBam(string, int, unsigned int) nogil except +
        int find_junctions(int min_experiments) nogil
        int find_junctions_from_region(vector[Gene*] gobj) nogil
        int ParseJunctionsFromFile() nogil
        int boostrap_samples(int msamples, int ksamples, np.float32_t * boots) nogil
        int get_njuncs() nogil
        map[string, unsigned int] get_junc_map() nogil

