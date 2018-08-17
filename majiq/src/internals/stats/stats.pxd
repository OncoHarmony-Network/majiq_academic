from libcpp.vector cimport vector

cdef extern from "stats.hpp":
    cdef cppclass HetStats:
        cdef bool initialize_statistics(vector<string> list_stats) nogil ;