from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "stats/stats.hpp":
    cdef cppclass HetStats:
        bint initialize_statistics(vector[string] list_stats) nogil ;