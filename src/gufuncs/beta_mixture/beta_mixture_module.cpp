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

#include "GlobalRNGPool.hpp"

#include "Approximation.hpp"
#include "CDF.hpp"
#include "DPsiDiscrete.hpp"
#include "Moments.hpp"
#include "PDF.hpp"
#include "PMF.hpp"
#include "Quantile.hpp"
#include "Sample.hpp"
#include "TTestSample.hpp"
#include "StatsSample.hpp"

// for functions requiring sampling
static PyObject* SetSeedGlobalGen(PyObject* self, PyObject* args) {
  int64_t seed;
  if (!PyArg_ParseTuple(args, "L", &seed)) {
    // failed to parse arguments, should raise exception in Python
    return nullptr;
  }
  MajiqGufuncs::BetaMixture::global_rng_pool.seed(seed);
  Py_RETURN_NONE;
}
static char SetSeedGlobalGen_doc[] = R"pbdoc(
rng_seed(seed: int) -> None

Set seed for pool of random number generators

Parameters
----------
seed: int
    seed random number generator with this value
)pbdoc";
static PyObject* IncreasePoolSizeGlobalGen(PyObject* self, PyObject* args) {
  int64_t n;
  if (!PyArg_ParseTuple(args, "L", &n)) {
    // failed to parse arguments, should raise exception in Python
    return nullptr;
  }
  MajiqGufuncs::BetaMixture::global_rng_pool.resize(n);
  Py_RETURN_NONE;
}
static char IncreasePoolSizeGlobalGen_doc[] = R"pbdoc(
rng_resize(n: int) -> None

Resize pool of random number generators to allow that many simultaneous threads

Parameters
----------
n: int
    minimum number of random number generators there should be available to work
    with
)pbdoc";

// identify available statistics
static PyObject* AvailableStatistics(PyObject* self, PyObject* args) {
  // no arguments
  if (!PyArg_ParseTuple(args, "")) {
    return nullptr;
  }
  // construct dictionary to return
  using MajiqGufuncs::BetaMixture::StatsSample::HetStats;
  return Py_BuildValue(
      "{s:L,s:L,s:L,s:L}",
      "ttest", static_cast<int64_t>(HetStats::TTest),
      "mannwhitneyu", static_cast<int64_t>(HetStats::MannWhitney),
      "tnom", static_cast<int64_t>(HetStats::TNOM),
      "infoscore", static_cast<int64_t>(HetStats::InfoScore));
}
static char AvailableStatistics_doc[] = R"pbdoc(
stats_available() -> Dict[str, int]

Mapping of available statistics to integer-value for use in stats_sample
)pbdoc";

// for functions with no extra data being passed in
static void *data[1] = {NULL};

// define module
static PyMethodDef ModuleMethods[] = {
  {"rng_seed", SetSeedGlobalGen, METH_VARARGS, SetSeedGlobalGen_doc},
  {"rng_resize", IncreasePoolSizeGlobalGen, METH_VARARGS,
    IncreasePoolSizeGlobalGen_doc},
  {"stats_available", AvailableStatistics, METH_VARARGS,
    AvailableStatistics_doc},
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

  namespace Approximation = MajiqGufuncs::BetaMixture::Approximation;
  PyObject *approximation = PyUFunc_FromFuncAndDataAndSignature(
      Approximation::funcs, data, Approximation::types,
      Approximation::ntypes, Approximation::nin, Approximation::nout,
      PyUFunc_None, Approximation::name, Approximation::doc, 0,
      Approximation::signature);
  PyDict_SetItemString(d, Approximation::name, approximation);
  Py_DECREF(approximation);

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

  namespace Sample = MajiqGufuncs::BetaMixture::Sample;
  PyObject *sample = PyUFunc_FromFuncAndDataAndSignature(
      Sample::funcs, data, Sample::types,
      Sample::ntypes, Sample::nin, Sample::nout,
      PyUFunc_None, Sample::name, Sample::doc, 0,
      Sample::signature);
  PyDict_SetItemString(d, Sample::name, sample);
  Py_DECREF(sample);

  namespace StatsSample = MajiqGufuncs::BetaMixture::StatsSample;
  PyObject *stats_sample = PyUFunc_FromFuncAndDataAndSignature(
      StatsSample::funcs, data, StatsSample::types,
      StatsSample::ntypes, StatsSample::nin, StatsSample::nout,
      PyUFunc_None, StatsSample::name, StatsSample::doc, 0,
      StatsSample::signature);
  PyDict_SetItemString(d, StatsSample::name, stats_sample);
  Py_DECREF(stats_sample);

  namespace TTestSample = MajiqGufuncs::BetaMixture::TTestSample;
  PyObject *ttest_sample = PyUFunc_FromFuncAndDataAndSignature(
      TTestSample::funcs, data, TTestSample::types,
      TTestSample::ntypes, TTestSample::nin, TTestSample::nout,
      PyUFunc_None, TTestSample::name, TTestSample::doc, 0,
      TTestSample::signature);
  PyDict_SetItemString(d, TTestSample::name, ttest_sample);
  Py_DECREF(ttest_sample);

  namespace Quantile = MajiqGufuncs::BetaMixture::Quantile;
  PyObject *quantile = PyUFunc_FromFuncAndDataAndSignature(
      Quantile::funcs, data, Quantile::types,
      Quantile::ntypes, Quantile::nin, Quantile::nout,
      PyUFunc_None, Quantile::name, Quantile::doc, 0,
      Quantile::signature);
  PyDict_SetItemString(d, Quantile::name, quantile);
  Py_DECREF(quantile);

  namespace DPsiDiscrete = MajiqGufuncs::BetaMixture::DPsiDiscrete;
  PyObject *dpsi_discrete = PyUFunc_FromFuncAndDataAndSignature(
      DPsiDiscrete::funcs, data, DPsiDiscrete::types,
      DPsiDiscrete::ntypes, DPsiDiscrete::nin, DPsiDiscrete::nout,
      PyUFunc_None, DPsiDiscrete::name, DPsiDiscrete::doc, 0,
      DPsiDiscrete::signature);
  PyDict_SetItemString(d, DPsiDiscrete::name, dpsi_discrete);
  Py_DECREF(dpsi_discrete);

  // return pointer to final module object
  return m;
}
