/**
 * Gufunc implementation of nanmedian in MAJIQ
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Authors: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_NANMEDIAN_HPP
#define MAJIQGUFUNCS_NANMEDIAN_HPP

#include <numpy/ndarraytypes.h>

#include <limits>
#include <vector>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/quantile.hpp>

namespace MajiqGufuncs {
namespace NanMedian {

static char name[] = "nanmedian";
constexpr int nin = 1;
constexpr int nout = 1;
static char signature[] = "(n)->()";
static char doc[] = R"pbdoc(
Compute medians omitting nans

Parameters
----------
x: array[float]
  Input array with values to take median over on core axis, ignoring nans

Returns
-------
array[float]
  medians over core axis
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
  const npy_intp dim_x = dimensions[0];
  // inner strides
  const npy_intp str_x = steps[0];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto x = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto out = CoreIt<RealT>::begin(args[1], outer_stride[1]);

  // trivial cases
  if (dim_broadcast < 1) { return; }
  if (dim_x < 1) {
    // no values of x, so nan
    out.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    return;
  }
  if (dim_x == 1) {
    // handle case with dim_broadcast <= 1 because zero stride --> infinite loop
    if (dim_broadcast > 1) {
      std::copy(x, x + dim_broadcast, out);
    } else {
      *out = *x;
    }
    return;
  }
  // otherwise, more than one value of x
  // buffer for values of x that we can mutate
  std::vector<RealT> x_notnan(dim_x);
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++out) {
    auto x_n = x.with_stride(str_x);
    // copy over not-nan values of x
    auto x_notnan_end = std::copy_if(x_n, x_n + dim_x, x_notnan.begin(),
        [](RealT y) { return !std::isnan(y); });
    using MajiqInclude::median;
    *out = median(x_notnan.begin(), x_notnan_end);
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
  NPY_FLOAT, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace NanMedian
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_NANMEDIAN_HPP
