from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from majiq.src.constants import *
from majiq.src.internals.HetStats cimport HetStats
from majiq.src.internals.qLSV cimport dpsiLSV, hetLSV, qLSV, psiLSV
from majiq.src.internals.mtypes cimport *
import numpy as np
import pickle
import sys
cimport numpy as np
import cython

cdef extern from "psi.hpp":

    cdef psi_distr_t& get_psi_border(psi_distr_t& psi_border, int nbins) nogil ;
    # cdef psi_distr_t& get_psi_border(int nbins) nogil ;
    cdef void psi_posterior(psiLSV* lsvObj, psi_distr_t& psi_border, int nbins) nogil ;

    cdef void deltapsi_posterior(dpsiLSV* lsvObj, vector[psi_distr_t]& prior_matrix, psi_distr_t& psi_border,
                                 int nbins) nogil ;

    cdef void get_samples_from_psi2(vector[psi_distr_t]& i_psi, np.float32_t* osamps, np.float32_t* o_mu_psi,
                                   np.float32_t* o_postpsi, int psi_samples, int j_offset, psi_distr_t& psi_border, int njunc,
                                   int msamples, int nbins, bint is_ir) nogil ;

    cdef void get_samples_from_psi(float* osamps, hetLSV* lsvObj, int psi_samples, psi_distr_t psi_border,
                                   int nbins, int cidx, int fidx) nogil ;

    cdef void get_samples_from_psi3(vector[psi_distr_t]& i_psi, vector[psi_distr_t]& osamps, psi_distr_t& o_mupsi,
                                   vector[psi_distr_t]& o_postpsi, int psi_samples, int j_offset,
                                   psi_distr_t psi_border, int njunc, int msamples, int nbins, bint is_ir) nogil ;

    cdef void test_calc(vector[psi_distr_t]& oPvals, HetStats* HetStatsObj, hetLSV* lsvObj, int psamples,
                        np.float32_t quant) nogil ;

    cdef void adjustdelta(psi_distr_t& o_mixtpdf, psi_distr_t& emp_dpsi, int num_iter, int nbins) nogil ;



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
    # return delta_psi, delta_psi_ir

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef inline void __load_default_prior(vector[vector[psi_distr_t]]& prior_matrix):

    cdef int numbins = prior_matrix[0].size()
    cdef int xx, yy

    encoding = sys.getfilesystemencoding()
    direc = os.path.dirname(__file__)

    fop = open('%s/../data/defaultprior.pickle' % direc, 'rb')
    fast_pickler = pickle.Unpickler(fop)
    data = fast_pickler.load().astype(np.float32)
    fop.close()
    # print_prior(data, numbins)
    data /= np.sum(data)
    # for xx in range(numbins):
    #     for yy in range(numbins):
    #         prior_matrix[0][xx][yy] = np.log(data[xx][yy])
    #         prior_matrix[1][xx][yy] = np.log(data[xx][yy])

    return

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef inline print_prior(vector[vector[psi_distr_t]] matrix, int nbins):
    cdef float sum = 0.0
    for xx in range(nbins):
        for yy in range(nbins):
            sys.stderr.write("%.4f, " % np.exp(matrix[0][xx][yy]))
            sum += np.exp(matrix[0][xx][yy])
        sys.stderr.write("\n")
    sys.stderr.write('##MATRIX [0] sum: %.4f\n' % sum)

    sum = 0
    for xx in range(nbins):
        for yy in range(nbins):
            sys.stderr.write("%.4f, " % np.exp(matrix[1][xx][yy]))
            sum += np.exp(matrix[1][xx][yy])
        sys.stderr.write("\n")
    sys.stderr.write('##MATRIX [1] sum: %.4f\n' % sum)


cdef inline void gen_prior_matrix(vector[vector[psi_distr_t]]& prior_matrix, dict lsv_type, dict lsv_empirical_psi1,
                            dict lsv_empirical_psi2, str output, list names, str plotpath, int iter, float binsize,
                            int numbins, bint defaultprior, int minpercent, object logger):

    cdef psi_distr_t mixture_pdf = psi_distr_t(numbins*2)
    cdef list list_of_lsv, njun_prior
    cdef int prior_idx
    cdef np.ndarray[np.float32_t, ndim=1] best_deltap, best_dpsi, best_dpsi_ir
    cdef np.ndarray[np.float32_t, ndim=1] best_delta_psi
    cdef np.ndarray[np.float32_t, ndim=3] np_pmatrix = np.zeros(shape=(2, numbins, numbins), dtype=np.float32)


    #Start prior matrix
    logger.info("Calculating prior matrix...")
    if defaultprior:
        __load_default_prior(prior_matrix)
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
                    __load_default_prior(prior_matrix)
                else:
                    prior_matrix[1] = prior_matrix[0]
                break

            logger.debug("Parametrizing 'best set'...%s", prior_idx)
            adjustdelta(mixture_pdf, best_delta_psi, iter, numbins*2)
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

        # plot_matrix(prior_matrix[prior_idx], "Prior Matrix , version %s" % prior_idx,
        #             "prior_matrix_jun_%s" % nj, plotpath)

    for xx in range(numbins):
        for yy in range(numbins):
            prior_matrix[0][xx][yy] = np.log(np_pmatrix[0, xx, yy])
            prior_matrix[1][xx][yy] = np.log(np_pmatrix[1, xx, yy])
            # print('KLKKK2')
    # print_prior(prior_ma