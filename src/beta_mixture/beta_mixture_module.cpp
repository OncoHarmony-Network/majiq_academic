/**
 * beta_mixture_module.cpp
 *
 * Defines Python module with Numpy gufuncs for mixture of beta distributions
 *
 * Author: Joseph K Aicher
 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/npy_3kcompat.h"
#include "cdf.h"
#include "pdf.h"
#include "quantile.h"
#include "moments.h"

using namespace BetaMixture;

// all functions have the same types
static char types[8] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE
};


// no extra data being passed in
static void *data[1] = {NULL};


// functions
PyUFuncGenericFunction funcs_moments[2] = {
    (PyUFuncGenericFunction) &(mixture_moments_outer<npy_float>),
    (PyUFuncGenericFunction) &(mixture_moments_outer<npy_double>)
};
const char *signature_moments = "(n),(n)->(),()";
const char *name_moments = "moments";
const char *doc_moments = R"mtdoc(
        Mean and variance of uniform mixture of beta distributions

        Usage: mean, variance = moments(a, b)
        Signature: (n),(n)->(),()

        Computes mean and variance using parameters along last axis of a and b
        using numpy broadcasting.

        Parameters
        ----------
        x1, x2: array_like
            arrays of beta distribution parameters. Must share last axis,
            otherwise normal numpy broadcasting rules apply

        Returns
        -------
        mean, variance: Tuple[ndarray, ndarray]
            Tuple with mean and variance of mixture distributions
        )mtdoc";

PyUFuncGenericFunction funcs_cdf[2] = {
    (PyUFuncGenericFunction) &(mixture_cdf_outer<npy_float>),
    (PyUFuncGenericFunction) &(mixture_cdf_outer<npy_double>)
};
const char *signature_cdf = "(),(n),(n)->()";
const char *name_cdf = "cdf";
const char *doc_cdf = R"mtdoc(
        Cumulative distribution function of uniform mixture of beta distributions

        Usage: cdf(x, a, b)
        Signature: (),(n),(n)->()

        Computes values of CDF with parameters along last axis of a and b using
        numpy broadcasting.

        Parameters
        ----------
        x1: array_like
            locations in domain of distribution to evaluate CDF for broadcasted
            distribution parameters
        x2, x3: array_like
            arrays of beta distribution parameters. Must share last axis,
            otherwise normal numpy broadcasting rules apply

        Returns
        -------
        cdf: ndarray
            Probability for CDF of mixture distribution at evaluated points for
            each mixture distribution
        )mtdoc";

PyUFuncGenericFunction funcs_pdf[2] = {
    (PyUFuncGenericFunction) &(mixture_pdf_outer<npy_float>),
    (PyUFuncGenericFunction) &(mixture_pdf_outer<npy_double>)
};
const char *signature_pdf = "(),(n),(n)->()";
const char *name_pdf = "pdf";
const char *doc_pdf = R"mtdoc(
        Probability density function of uniform mixture of beta distributions

        Usage: pdf(x, a, b)
        Signature: (),(n),(n)->()

        Computes values of PDF with parameters along last axis of a and b using
        numpy broadcasting.

        Parameters
        ----------
        x1: array_like
            locations in domain of distribution to evaluate PDF for broadcasted
            distribution parameters
        x2, x3: array_like
            arrays of beta distribution parameters. Must share last axis,
            otherwise normal numpy broadcasting rules apply

        Returns
        -------
        pdf: ndarray
            Probability density of mixture distribution at evaluated points for
            each mixture distribution
        )mtdoc";

PyUFuncGenericFunction funcs_quantile[2] = {
    (PyUFuncGenericFunction) &(mixture_quantile_outer<npy_float>),
    (PyUFuncGenericFunction) &(mixture_quantile_outer<npy_double>)
};
const char *signature_quantile = "(),(n),(n)->()";
const char *name_quantile = "quantile";
const char *doc_quantile = R"mtdoc(
        Quantile function of uniform mixture of beta distributions

        Usage: quantile(q, a, b)
        Signature: (),(n),(n)->()

        Computes quantiles of beta mixture distribution with parameters along
        last axis of a and b using numpy broadcasting.

        Parameters
        ----------
        x1: array_like
            Probabilities for which quantiles will be evaluated in the domain
            of beta distribution mixtures for broadcasted distribution
            parameters
        x2, x3: array_like
            arrays of beta distribution parameters. Must share last axis,
            otherwise normal numpy broadcasting rules apply

        Returns
        -------
        quantile: ndarray
            Quantiles of mixture distribution at evaluated points for each
            mixture distribution

        Notes
        -----
        Quantiles are estimated by performing Newton iteration on a bracketed
        interval guaranteed to have the true result by Chebyshev's inequality
        )mtdoc";


// define module
static PyMethodDef ModuleMethods[] = {
    // gufuncs have to be initizlied elsewhere
    {NULL, NULL, 0, NULL}
};
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    // namespace for module
    "beta_mixture",
    NULL,
    -1,
    ModuleMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_beta_mixture(void) {
    PyObject *m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    PyObject *d = PyModule_GetDict(m);

    import_array();
    import_umath();

    PyObject *mixture_moments = PyUFunc_FromFuncAndDataAndSignature(
            funcs_moments, data, types, 2, 2, 2, PyUFunc_None,
            name_moments, doc_moments, 0, signature_moments);
    PyDict_SetItemString(d, name_moments, mixture_moments);
    Py_DECREF(mixture_moments);

    PyObject *mixture_cdf = PyUFunc_FromFuncAndDataAndSignature(
            funcs_cdf, data, types, 2, 3, 1, PyUFunc_None,
            name_cdf, doc_cdf, 0, signature_cdf);
    PyDict_SetItemString(d, name_cdf, mixture_cdf);
    Py_DECREF(mixture_cdf);

    PyObject *mixture_pdf = PyUFunc_FromFuncAndDataAndSignature(
            funcs_pdf, data, types, 2, 3, 1, PyUFunc_None,
            name_pdf, doc_pdf, 0, signature_pdf);
    PyDict_SetItemString(d, name_pdf, mixture_pdf);
    Py_DECREF(mixture_pdf);

    PyObject *mixture_quantile = PyUFunc_FromFuncAndDataAndSignature(
            funcs_quantile, data, types, 2, 3, 1, PyUFunc_None,
            name_quantile, doc_quantile, 0, signature_quantile);
    PyDict_SetItemString(d, name_quantile, mixture_quantile);
    Py_DECREF(mixture_quantile);

    return m;
}
