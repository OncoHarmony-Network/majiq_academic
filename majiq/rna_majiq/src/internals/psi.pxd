from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from rna_majiq.src.constants import *
from rna_majiq.src.internals.HetStats cimport HetStats
from rna_majiq.src.internals.qLSV cimport dpsiLSV, hetLSV, qLSV, psiLSV
from rna_majiq.src.internals.mtypes cimport *
import numpy as np
import pickle
import sys
cimport numpy as np
import cython

cdef extern from "<random>" namespace "std":
    cdef cppclass mt19937:
        # example from <https://stackoverflow.com/questions/40976880/canonical-way-to-generate-random-numbers-in-cython>
        mt19937() # we need to define this constructor to stack allocate classes in Cython
        mt19937(unsigned int seed) # not worrying about matching the exact int type for seed
        void seed(unsigned int seed)

cdef extern from "psi.hpp":

    cdef psi_distr_t& get_psi_border(psi_distr_t& psi_border, int nbins) nogil ;
    # cdef psi_distr_t& get_psi_border(int nbins) nogil ;
    cdef void psi_posterior(psiLSV* lsvObj, psi_distr_t& psi_border, int nbins) nogil ;

    cdef void deltapsi_posterior(dpsiLSV* lsvObj, vector[psi_distr_t]& prior_matrix, psi_distr_t& psi_border,
                                 int nbins) nogil ;

    cdef void get_samples_from_psi(
        float* osamps, hetLSV* lsvObj, int psi_samples, int visualization_samples,
        psi_distr_t psi_border, int nbins, int cidx, int fidx, mt19937 &generator
    ) nogil;

    cdef void test_calc(vector[psi_distr_t]& mean_pvalues, vector[psi_distr_t]& sample_pvalues, psi_distr_t& oScore, HetStats* HetStatsObj, hetLSV* lsvObj,
                        int psamples, np.float32_t quant) nogil ;

    cdef int adjustdelta(psi_distr_t& o_mixtpdf, psi_distr_t& emp_dpsi, int num_iter, int nbins) nogil ;


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef inline tuple _empirical_delta_psi(list list_of_lsv, dict lsv_empirical_psi1, dict lsv_empirical_psi2, object lsv_type):
    """
    Simple PSI calculation without involving a dirichlet prior, coming from reads from junctions
    """
    cdef string lsv
    cdef list delta_psi = []
    cdef list delta_psi_ir = []

    for lsv in list_of_lsv:
        if lsv_type[lsv][1] > 2 : continue
        # Assuming that the type is the same in all the replicas and groups
        if lsv_type[lsv][0].endswith(b'i'):
            delta_psi_res = delta_psi_ir
        else:
            delta_psi_res = delta_psi
        delta_psi_res.append(lsv_empirical_psi1[lsv][0] - lsv_empirical_psi2[lsv][0])
        delta_psi_res.append(lsv_empirical_psi2[lsv][0] - lsv_empirical_psi1[lsv][0])

    return np.array(delta_psi, dtype=np.float32), np.array(delta_psi_ir, dtype=np.float32)


cdef inline int gen_prior_matrix(vector[vector[psi_distr_t]]& prior_matrix, dict lsv_type, dict lsv_empirical_psi1,
                            dict lsv_empirical_psi2, str output, list names, int iter, float binsize,
                            int numbins, bint defaultprior, int minpercent, object logger) except -1:

    import pickle

    cdef psi_distr_t mixture_pdf = psi_distr_t(numbins*2)
    cdef list list_of_lsv, njun_prior
    cdef int prior_idx
    cdef np.ndarray[np.float32_t, ndim=1] best_deltap, best_dpsi, best_dpsi_ir
    cdef np.ndarray[np.float32_t, ndim=1] best_delta_psi
    cdef np.ndarray[np.float32_t, ndim=3] np_pmatrix = np.zeros(shape=(2, numbins, numbins), dtype=np.float32)

    #Start prior matrix
    logger.info("Calculating prior matrix...")
    if defaultprior:

        # TODO --- repeated code, refactor into cython function
        with open('%s/../data/defaultprior.pickle' % os.path.dirname(__file__), 'rb') as f:
            fast_pickler = pickle.Unpickler(f)
            data = fast_pickler.load().astype(np.float32)

        data /= np.sum(data)
        for xx in range(numbins):
            for yy in range(numbins):
                np_pmatrix[0][xx][yy] = data[xx][yy]
                np_pmatrix[1][xx][yy] = data[xx][yy]
        # end repeated code

    else:

        logger.debug('Filtering to obtain "best set"...')

        list_of_lsv = list(set(lsv_empirical_psi1.keys()).intersection(set(lsv_empirical_psi2.keys())))
        logger.debug("'Best set' is %s events" % len(list_of_lsv))
        best_dpsi, best_dpsi_ir = _empirical_delta_psi(list_of_lsv, lsv_empirical_psi1, lsv_empirical_psi2, lsv_type)

        for prior_idx, best_delta_psi in enumerate((best_dpsi, best_dpsi_ir)):
            # initialize mixture_pdf
            for i in range(numbins*2):
                mixture_pdf[i] = 0

            if len(best_delta_psi) <= 100:
                if prior_idx == 0:

                    # TODO --- repeated code, refactor into cython function
                    with open('%s/../data/defaultprior.pickle' % os.path.dirname(__file__), 'rb') as f:
                        fast_pickler = pickle.Unpickler(f)
                        data = fast_pickler.load().astype(np.float32)

                    data /= np.sum(data)
                    for xx in range(numbins):
                        for yy in range(numbins):
                            np_pmatrix[0][xx][yy] = data[xx][yy]
                            np_pmatrix[1][xx][yy] = data[xx][yy]
                    # end repeated code

                else:
                    np_pmatrix[1][:] = np_pmatrix[0][:]
                break

            logger.debug("Parametrizing 'best set'...%s", prior_idx)
            r = adjustdelta(mixture_pdf, best_delta_psi, iter, numbins*2)
            if r == -1 :
                raise ValueError(" The input data does not have enought statistic power in order to calculate "
                                 "the prior. Check if the input is correct or use the --default-prior option in "
                                  " order to use a precomputed prior")
            for i in range(numbins):
                for j in range(numbins):
                    np_pmatrix[prior_idx][i][j] = mixture_pdf[j-i+(numbins-1)]

            if np.isnan(np_pmatrix[prior_idx]).any():
                if prior_idx == 1:
                    logger.warning("Not enought statistic power to calculate the intron retention specific prior, "
                                   "in that case we will use the global prior")
                    np_pmatrix[prior_idx] = np_pmatrix[0]
                else:
                    raise ValueError(" The input data does not have enought statistic power in order to calculate "
                                     "the prior. Check if the input is correct or use the --default-prior option in "
                                     " order to use a precomputed prior")
            else:
                np_pmatrix[prior_idx] /= np.sum(np_pmatrix[prior_idx])

            # renormalize so it sums 1

    for xx in range(numbins):
        for yy in range(numbins):
            prior_matrix[0][xx][yy] = np.log(np_pmatrix[0, xx, yy])
            prior_matrix[1][xx][yy] = np.log(np_pmatrix[1, xx, yy])
