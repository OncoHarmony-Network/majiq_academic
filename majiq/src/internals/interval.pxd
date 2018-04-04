cdef extern from "interval.hpp" namespace "interval":
    cdef struct Interval:
        int low
        int high
        void * ptr

    cdef struct ITNode:
        Interval *i;
        int max
        ITNode *left
        ITNode *right

    Interval *intervalSearch(ITNode *root, Interval i) nogil
    void inorder(ITNode *root) nogil
    ITNode *insert(ITNode *root, Interval i) nogil