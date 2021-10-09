/*
 * Implementation of Welch's t-test (two-sample t-test assuming unequal
 * variances) for use in MAJIQ
 *
 * Copyright 2020 <University of Pennsylvania>
 *
 * Authors: Joseph K Aicher, Jorge Vaquero-Garcia
 */

#ifndef MAJIQGUFUNCS_TTEST_HPP
#define MAJIQGUFUNCS_TTEST_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/stats/TTest.hpp>


namespace MajiqGufuncs {
namespace TTest {

static char name[] = "ttest";
constexpr int nin = 2;
constexpr int nout = 1;
static char signature[] = "(n),(n)->()";
static char doc[] = R"pbdoc(
Compute p-values for Welch's t-test on input data

Compute p-values for Welch's t-test on input data, using a two-sided
alternative hypothesis and omitting nan values.

Parameters
----------
x: array[float]
    test over observations in last axis
labels: array[bool]
    test over labels in last axis

Returns
-------
array[float]
    broadcast p-values for observations/labels. Invalid tests are nan.
)pbdoc";

template <typename RealT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides on each variable for outer loop
  const npy_intp* outer_stride = steps;
  steps += nin + nout;

  // core dimensions
  const npy_intp dim_core = dimensions[0];
  // inner strides
  const npy_intp inner_stride_x = steps[0];
  const npy_intp inner_stride_labels = steps[1];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto x = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto labels = CoreIt<bool>::begin(args[1], outer_stride[1]);
  auto out = CoreIt<RealT>::begin(args[2], outer_stride[2]);

  if (dim_core < 1) {
    // no samples to process, so must be nan
    out.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    return;
  }

  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++labels, ++out) {
    using MajiqInclude::TTest::Test;
    *out = Test(
        x.with_stride(inner_stride_x),
        labels.with_stride(inner_stride_labels),
        dim_core);
  }
  return;
}

constexpr int ntypes = 2;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>)
};
static char types[ntypes * (nin + nout)] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_BOOL, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_BOOL, NPY_DOUBLE,
};


}  // namespace TTest
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_TTEST_HPP
