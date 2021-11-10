/**
 * stats_module.cpp
 *
 * Helper gufuncs for working with statistics
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

#include "InfoScore.hpp"
#include "MannWhitney.hpp"
#include "TNOM.hpp"
#include "TTest.hpp"

#include "NanMedian.hpp"
#include "NanQuantile.hpp"


// no extra data being passed in
static void *data[1] = {NULL};


// define module
static PyMethodDef ModuleMethods[] = {
  {NULL, NULL, 0, NULL}
};
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  // namespace for module
  "_stats",
  NULL,
  -1,
  ModuleMethods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit__stats(void) {
  PyObject* m = PyModule_Create(&moduledef);
  if (!m) {
    return NULL;
  }
  // get instance dictionary for module, which we will attach gufuncs to
  PyObject* d = PyModule_GetDict(m);

  import_array();
  import_umath();

  /**
   * Perform statistical tests
   */

  namespace TTest = MajiqGufuncs::TTest;
  PyObject *ttest = PyUFunc_FromFuncAndDataAndSignature(
      TTest::funcs, data, TTest::types,
      TTest::ntypes, TTest::nin, TTest::nout,
      PyUFunc_None, TTest::name, TTest::doc, 0,
      TTest::signature);
  PyDict_SetItemString(d, TTest::name, ttest);
  Py_DECREF(ttest);

  namespace MannWhitney = MajiqGufuncs::MannWhitney;
  PyObject *mannwhitneyu = PyUFunc_FromFuncAndDataAndSignature(
      MannWhitney::funcs, data, MannWhitney::types,
      MannWhitney::ntypes, MannWhitney::nin, MannWhitney::nout,
      PyUFunc_None, MannWhitney::name, MannWhitney::doc, 0,
      MannWhitney::signature);
  PyDict_SetItemString(d, MannWhitney::name, mannwhitneyu);
  Py_DECREF(mannwhitneyu);

  namespace TNOM = MajiqGufuncs::TNOM;
  PyObject *tnom = PyUFunc_FromFuncAndDataAndSignature(
      TNOM::funcs, data, TNOM::types,
      TNOM::ntypes, TNOM::nin, TNOM::nout,
      PyUFunc_None, TNOM::name, TNOM::doc, 0,
      TNOM::signature);
  PyDict_SetItemString(d, TNOM::name, tnom);
  Py_DECREF(tnom);

  namespace InfoScore = MajiqGufuncs::InfoScore;
  PyObject *infoscore = PyUFunc_FromFuncAndDataAndSignature(
      InfoScore::funcs, data, InfoScore::types,
      InfoScore::ntypes, InfoScore::nin, InfoScore::nout,
      PyUFunc_None, InfoScore::name, InfoScore::doc, 0,
      InfoScore::signature);
  PyDict_SetItemString(d, InfoScore::name, infoscore);
  Py_DECREF(infoscore);

  namespace NanMedian = MajiqGufuncs::NanMedian;
  PyObject *nanmedian = PyUFunc_FromFuncAndDataAndSignature(
      NanMedian::funcs, data, NanMedian::types,
      NanMedian::ntypes, NanMedian::nin, NanMedian::nout,
      PyUFunc_None, NanMedian::name, NanMedian::doc, 0,
      NanMedian::signature);
  PyDict_SetItemString(d, NanMedian::name, nanmedian);
  Py_DECREF(nanmedian);

  namespace NanQuantile = MajiqGufuncs::NanQuantile;
  PyObject *nanquantile = PyUFunc_FromFuncAndDataAndSignature(
      NanQuantile::funcs, data, NanQuantile::types,
      NanQuantile::ntypes, NanQuantile::nin, NanQuantile::nout,
      PyUFunc_None, NanQuantile::name, NanQuantile::doc, 0,
      NanQuantile::signature);
  PyDict_SetItemString(d, NanQuantile::name, nanquantile);
  Py_DECREF(nanquantile);

  // return pointer to final module object
  return m;
}
