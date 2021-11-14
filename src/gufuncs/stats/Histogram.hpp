/**
 * Gufunc implementation of histogram over uniform bins in MAJIQ
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Authors: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_HISTOGRAM_HPP
#define MAJIQGUFUNCS_HISTOGRAM_HPP

#include <numpy/ndarraytypes.h>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/histogram.hpp>

namespace MajiqGufuncs {
namespace Histogram {

static char name[] = "histogram";
constexpr int nin = 4;
constexpr int nout = 1;
static char signature[] = "(n),(),(),(b)->(b)";
static char doc[] = R"pbdoc(
Compute histogram over uniform bins ignoring nans

Parameters
----------
x: array[float]
  Values to be summarized as histogram. NaN values are ignored
x_min, x_max: float
  Minimum and maximum values to be counted (values outside of range ignored)
dummy: array[float]
  Array with shape equal to desired number of output bins

Returns
-------
array[int]: counts of input values on uniform bins within [x_min, x_max)
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
  const npy_intp dim_out = dimensions[1];
  // inner strides
  const npy_intp str_x = steps[0];
  // ignore str_dummy = steps[1];
  const npy_intp str_out = steps[2];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto x = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto xmin = CoreIt<RealT>::begin(args[1], outer_stride[1]);
  auto xmax = CoreIt<RealT>::begin(args[2], outer_stride[2]);
  // ignoring dummy {args[3], outer_stride[3]}
  auto out = CoreIt<int64_t>::begin(args[4], outer_stride[4]);

  // trivial cases
  if (dim_out < 1) { return; }
  if (dim_x < 1) {
    // no values of x, so all counts are zero
    if (dim_out > 1) {
      for (npy_intp i = 0; i < dim_broadcast; ++i, ++out) {
        auto out_out = out.with_stride(str_out);
        std::fill(out_out, out_out + dim_out, int64_t{0});
      }
    } else {
      if (dim_broadcast > 1) {
        std::fill(out, out + dim_broadcast, int64_t{0});
      } else {
        *out = int64_t{0};
      }
    }
    return;
  }

  // outer loop over broadcast dimension
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++xmin, ++xmax, ++out) {
    using MajiqInclude::Histogram::histogram;
    histogram(x.with_stride(str_x), dim_x, *xmin, *xmax,
        out.with_stride(str_out), dim_out);
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
  NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_DOUBLE, NPY_INT64,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_INT64,
};

}  // namespace Histogram
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_HISTOGRAM_HPP
