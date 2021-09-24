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

#include "ClipAndOffsetSum.hpp"
#include "ClipAndNormalize.hpp"
#include "OffsetOr.hpp"

#include "InfoScore.hpp"
#include "MannWhitney.hpp"
#include "TNOM.hpp"
#include "TTest.hpp"


// no extra data being passed in
static void *data[1] = {NULL};


// define module
static PyMethodDef ModuleMethods[] = {
  {NULL, NULL, 0, NULL}
};
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  // namespace for module
  "gufuncs",
  NULL,
  -1,
  ModuleMethods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit_gufuncs(void) {
  PyObject* m = PyModule_Create(&moduledef);
  if (!m) {
    return NULL;
  }
  // get instance dictionary for module, which we will attach gufuncs to
  PyObject* d = PyModule_GetDict(m);

  import_array();
  import_umath();

  /**
   *Quickly work with groups defined by offsets
   */

  namespace OffsetOr = MajiqGufuncs::OffsetOr;
  PyObject *offset_logical_or = PyUFunc_FromFuncAndDataAndSignature(
      OffsetOr::funcs, data, OffsetOr::types,
      OffsetOr::ntypes, OffsetOr::nin, OffsetOr::nout,
      PyUFunc_None, OffsetOr::name, OffsetOr::doc, 0,
      OffsetOr::signature);
  PyDict_SetItemString(d, OffsetOr::name, offset_logical_or);
  Py_DECREF(offset_logical_or);

  namespace ClipAndOffsetSum = MajiqGufuncs::ClipAndOffsetSum;
  PyObject *clip_and_offsetsum = PyUFunc_FromFuncAndDataAndSignature(
      ClipAndOffsetSum::funcs, data, ClipAndOffsetSum::types,
      ClipAndOffsetSum::ntypes, ClipAndOffsetSum::nin, ClipAndOffsetSum::nout,
      PyUFunc_None, ClipAndOffsetSum::name, ClipAndOffsetSum::doc, 0,
      ClipAndOffsetSum::signature);
  PyDict_SetItemString(d, ClipAndOffsetSum::name, clip_and_offsetsum);
  Py_DECREF(clip_and_offsetsum);

  PyObject *offsetsum = PyUFunc_FromFuncAndDataAndSignature(
      ClipAndOffsetSum::noclipfuncs, data, ClipAndOffsetSum::types,
      ClipAndOffsetSum::ntypes, ClipAndOffsetSum::nin, ClipAndOffsetSum::nout,
      PyUFunc_None, ClipAndOffsetSum::noclipname, ClipAndOffsetSum::noclipdoc,
      0, ClipAndOffsetSum::signature);
  PyDict_SetItemString(d, ClipAndOffsetSum::noclipname, offsetsum);
  Py_DECREF(offsetsum);

  namespace ClipAndNormalize = MajiqGufuncs::ClipAndNormalize;
  PyObject *clip_and_normalize = PyUFunc_FromFuncAndDataAndSignature(
      ClipAndNormalize::funcs, data, ClipAndNormalize::types,
      ClipAndNormalize::ntypes, ClipAndNormalize::nin, ClipAndNormalize::nout,
      PyUFunc_None, ClipAndNormalize::name, ClipAndNormalize::doc, 0,
      ClipAndNormalize::signature);
  PyDict_SetItemString(d, ClipAndNormalize::name, clip_and_normalize);
  Py_DECREF(clip_and_normalize);

  PyObject *clip_and_normalize_strict = PyUFunc_FromFuncAndDataAndSignature(
      ClipAndNormalize::strictfuncs, data, ClipAndNormalize::types,
      ClipAndNormalize::ntypes, ClipAndNormalize::nin, ClipAndNormalize::nout,
      PyUFunc_None, ClipAndNormalize::strictname, ClipAndNormalize::strictdoc,
      0, ClipAndNormalize::signature);
  PyDict_SetItemString(
      d, ClipAndNormalize::strictname, clip_and_normalize_strict);
  Py_DECREF(clip_and_normalize_strict);


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

  // return pointer to final module object
  return m;
}
