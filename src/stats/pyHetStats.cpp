/**
 * pyHetStats.cpp
 *
 * Set up test statistics for MAJIQ HET
 *
 * Copyright 2021 <University of Pennsylvania
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdexcept>


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;


/**
 * Compute p-values for Mann-Whitney U test on input data
 *
 * @param x: 2D array, test for each row over observations in columns
 * @param sortx: 2D array, indexes that sort x per row
 * @param labels: 2D array, class labels for each observation
 *
 * @return p-value for test statistic for each row (1D array)
 */
template <typename RealT>
py::array_t<RealT> MannWhitneyU(
    py::array_t<RealT> x,
    py::array_t<ssize_t> sortx,
    py::array_t<bool> labels) {
  // check for correct number of dimensions
  if (x.ndim() != 2) {
    throw std::runtime_error("x is not 2-dimensional");
  } else if (sortx.ndim() != 2) {
    throw std::runtime_error("sortx is not 2-dimensional");
  } else if (labels.ndim() != 2) {
    throw std::runtime_error("labels is not 2-dimensional");
  }
  // check that dimensions match
  if (x.shape(0) != sortx.shape(0) || x.shape(1) != sortx.shape(1)) {
    throw std::runtime_error("x.shape does not match sortx.shape");
  } else if (x.shape(0) != labels.shape(0) || x.shape(1) != labels.shape(1)) {
    throw std::runtime_error("x.shape does not match labels.shape");
  }

  // create output array, 1D with length equal to rows of x
  py::array_t<RealT> result(x.shape(0));
  // unchecked access to the array values
  auto _x = x.template unchecked<2>();
  auto _sortx = sortx.template unchecked<2>();
  auto _labels = labels.template unchecked<2>();
  auto _result = result.template mutable_unchecked<1>();

  // calculate statistics per row of x
  for (py::ssize_t i = 0; i < _x.shape(0); ++i) {
    // TODO(jaicher): actually do calculation
    _result(i) = RealT{1.};
  }

  // return result
  return result;
}


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
      &MannWhitneyU<float>, pbdoc_mannwhitneyu,
      py::arg("x"), py::arg("sortx"), py::arg("labels"));
  m.def("mannwhitneyu",
      &MannWhitneyU<double>, pbdoc_mannwhitneyu,
      py::arg("x"), py::arg("sortx"), py::arg("labels"));


  // end define functions

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}  // PYBIND11_MODULE
