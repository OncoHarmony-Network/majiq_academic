from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.set cimport set

cdef extern from "stats/stats.hpp":
    cdef cppclass HetStats:
        vector[string] names ;
        bint initialize_statistics(vector[string] list_stats) nogil ;
        int get_number_stats() nogil ;