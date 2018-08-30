from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.set cimport set

cdef extern from "stats/stats.hpp":
    cdef cppclass HetStats:
        set[string] names ;
        bint initialize_statistics(vector[string] list_stats) nogil ;