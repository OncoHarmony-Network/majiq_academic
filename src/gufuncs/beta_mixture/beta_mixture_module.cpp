/**
 * gufuncs_module.cpp
 *
 * Helper gufuncs for working with LSVs/PSI
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */


#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>

#include "CDF.hpp"
#include "Moments.hpp"
#include "PDF.hpp"
#include "PMF.hpp"
#include "Quantile.hpp"


// no extra data being passed in
static void *data[1] = {NULL};


// define module
static PyMethodDef ModuleMethods[] = {
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
  PyObject* m = PyModule_Create(&moduledef);
  if (!m) {
    return NULL;
  }
  // get instance dictionary for module, which we will attach gufuncs to
  PyObject* d = PyModule_GetDict(m);

  import_array();
  import_umath();

  namespace Moments = MajiqGufuncs::BetaMixture::Moments;
  PyObject *moments = PyUFunc_FromFuncAndDataAndSignature(
      Moments::funcs, data, Moments::types,
      Moments::ntypes, Moments::nin, Moments::nout,
      PyUFunc_None, Moments::name, Moments::doc, 0,
      Moments::signature);
  PyDict_SetItemString(d, Moments::name, moments);
  Py_DECREF(moments);

  namespace CDF = MajiqGufuncs::BetaMixture::CDF;
  PyObject *cdf = PyUFunc_FromFuncAndDataAndSignature(
      CDF::funcs, data, CDF::types,
      CDF::ntypes, CDF::nin, CDF::nout,
      PyUFunc_None, CDF::name, CDF::doc, 0,
      CDF::signature);
  PyDict_SetItemString(d, CDF::name, cdf);
  Py_DECREF(cdf);

  namespace PDF = MajiqGufuncs::BetaMixture::PDF;
  PyObject *pdf = PyUFunc_FromFuncAndDataAndSignature(
      PDF::funcs, data, PDF::types,
      PDF::ntypes, PDF::nin, PDF::nout,
      PyUFunc_None, PDF::name, PDF::doc, 0,
      PDF::signature);
  PyDict_SetItemString(d, PDF::name, pdf);
  Py_DECREF(pdf);

  namespace PMF = MajiqGufuncs::BetaMixture::PMF;
  PyObject *pmf = PyUFunc_FromFuncAndDataAndSignature(
      PMF::funcs, data, PMF::types,
      PMF::ntypes, PMF::nin, PMF::nout,
      PyUFunc_None, PMF::name, PMF::doc, 0,
      PMF::signature);
  PyDict_SetItemString(d, PMF::name, pmf);
  Py_DECREF(pmf);

  namespace Quantile = MajiqGufuncs::BetaMixture::Quantile;
  PyObject *quantile = PyUFunc_FromFuncAndDataAndSignature(
      Quantile::funcs, data, Quantile::types,
      Quantile::ntypes, Quantile::nin, Quantile::nout,
      PyUFunc_None, Quantile::name, Quantile::doc, 0,
      Quantile::signature);
  PyDict_SetItemString(d, Quantile::name, quantile);
  Py_DECREF(quantile);

  // return pointer to final module object
  return m;
}
