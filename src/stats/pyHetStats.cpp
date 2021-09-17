/**
 * pyHetStats.cpp
 *
 * Set up test statistics for MAJIQ HET
 *
 * Copyright 2021 <University of Pennsylvania
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "MannWhitney.hpp"
#include "TNOM.hpp"
#include "TTest.hpp"


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;




PYBIND11_MODULE(_stats, m) {
  m.attr("__name__") = "new_majiq._stats";
  // documentation of module
  m.doc() = R"pbdoc(
  MAJIQ implementations of two-sample statistical tests
  .. current_module::new_majiq._stats
  .. autosummary::
    ;toctree: _generate

  )pbdoc";

  // define functions

  // ttest
  constexpr char pbdoc_ttest[] = R"pbdoc(
  Compute p-values for Welch's t-test on input data

  Parameters
  ----------
  x: 2D array[float]
      test for each row over observations in columns
  labels: 2D array[bool]
      boolean class labels for each observation

  Returns
  -------
  1D array[float]
      p-values of test statistics for each row
  )pbdoc";
  m.def("ttest",
      &MajiqStats::TTest<float>, pbdoc_ttest,
      py::arg("x"), py::arg("labels"));
  m.def("ttest",
      &MajiqStats::TTest<double>, pbdoc_ttest,
      py::arg("x"), py::arg("labels"));

  // TNOM
  constexpr char pbdoc_tnom[] = R"pbdoc(
  Compute p-values for TNOM test on input data

  Parameters
  ----------
  x: 2D array[float]
      test for each row over observations in columns
  sortx: 2D array[int]
      indexes that sort x per row (i.e. np.argsort(x, axis=-1))
  labels: 2D array[bool]
      boolean class labels for each observation

  Returns
  -------
  1D array[float]
      p-values of test statistics for each row
  )pbdoc";
  m.def("tnom",
      &MajiqStats::TNOM<float>, pbdoc_tnom,
      py::arg("x"), py::arg("sortx"), py::arg("labels"));
  m.def("tnom",
      &MajiqStats::TNOM<double>, pbdoc_tnom,
      py::arg("x"), py::arg("sortx"), py::arg("labels"));

  // mannwhitneyu
  constexpr char pbdoc_mannwhitneyu[] = R"pbdoc(
  Compute p-values for Mann-Whitney U test on input data

  Parameters
  ----------
  x: 2D array[float]
      test for each row over observations in columns
  sortx: 2D array[int]
      indexes that sort x per row (i.e. np.argsort(x, axis=-1))
  labels: 2D array[bool]
      boolean class labels for each observation

  Returns
  -------
  1D array[float]
      p-values of test statistics for each row
  )pbdoc";
  m.def("mannwhitneyu",
      &MajiqStats::MannWhitneyU<float>, pbdoc_mannwhitneyu,
      py::arg("x"), py::arg("sortx"), py::arg("labels"));
  m.def("mannwhitneyu",
      &MajiqStats::MannWhitneyU<double>, pbdoc_mannwhitneyu,
      py::arg("x"), py::arg("sortx"), py::arg("labels"));


  // end define functions

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}  // PYBIND11_MODULE
