from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
cimport numpy as np
from majiq.src.internals.interval cimport ITNode
from majiq.src.internals.grimoire cimport Junction, Gene

ctypedef np.float64_t DTYPE_t
ctypedef vector[Gene *] gene_vect
cdef extern from "io_bam.hpp" namespace "io_bam":
    cdef cppclass IOBam:
        IOBam() nogil except +
        IOBam(string, int, unsigned int, unsigned int, map[string, gene_vect]) nogil except +
        IOBam(string, int, unsigned int) nogil except +
        # int find_junctions(int min_experiments) nogil
        # int find_junctions_from_region(vector[Gene*] gobj) nogil
        int ParseJunctionsFromFile() nogil
        # int ParseJunctionsFromFile(map[string, ITNode*]) nogil
        int boostrap_samples(int msamples, int ksamples, np.float32_t * boots) nogil
        int get_njuncs() nogil
        map[string, unsigned int] get_junc_map() nogil
        vector[Junction *] get_junc_vec() nogil
        int * get_junc_vec_summary() nogil ;

