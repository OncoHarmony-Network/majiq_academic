/**
 * Gufunc implementation of nanquantile (without interpolation options) for use
 * in MAJIQ (numpy implementation is not vectorized)
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Authors: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_NANQUANTILE_HPP
#define MAJIQGUFUNCS_NANQUANTILE_HPP

#include <numpy/ndarraytypes.h>

#include <limits>
#include <vector>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/quantile.hpp>

namespace MajiqGufuncs {
namespace NanQuantile {

static char name[] = "nanquantile";
constexpr int nin = 2;
constexpr int nout = 1;
static char signature[] = "(n),(q)->(q)";
static char doc[] = R"pbdoc(
Compute quantiles omitting nans (using linear interpolation)

Parameters
----------
x: array[float]
  Input array with values to take quantiles over on core axis, ignoring nans
q: array[float]
  Quantile or sequence of quantiles to compute, which must be between 0 and 1
  inclusive. NaN values will return NaN, values outside of range will be
  treated as 0 or 1 (whichever is closer).

Returns
-------
array[float]
  broadcast quantiles over core axis after ignoring nans
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
  const npy_intp dim_q = dimensions[1];
  // inner strides
  const npy_intp str_x = steps[0];
  const npy_intp str_q = steps[1];
  const npy_intp str_out = steps[2];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto x = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto q = CoreIt<npy_double>::begin(args[1], outer_stride[1]);
  auto out = CoreIt<RealT>::begin(args[2], outer_stride[2]);

  // handle trivial cases
  if (dim_q < 1) {
    return;  // out has zero dimension
  }
  if (dim_x < 1) {
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++out) {
      // no values of x, so nan
      out
        .with_stride(str_out)
        .fill(dim_q, std::numeric_limits<RealT>::quiet_NaN());
    }
    return;
  } else if (dim_x == 1) {
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++out) {
      // only one value of x, so whatever that value is
      out.with_stride(str_out).fill(dim_q, *x);
    }
    return;
  }
  // otherwise, more than one value of x and at least one requested quantile
  // buffer for values of x that we can mutate
  std::vector<RealT> x_notnan(dim_x);
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++q, ++out) {
    auto x_n = x.with_stride(str_x);
    auto q_q = q.with_stride(str_q);
    auto out_q = out.with_stride(str_out);
    // copy over not-nan values of x
    auto x_notnan_end = std::copy_if(x_n, x_n + dim_x, x_notnan.begin(),
        [](RealT y) { return !std::isnan(y); });
    // get quantiles
    for (npy_intp j = 0; j < dim_q; ++j, ++q_q, ++out_q) {
      if (std::isnan(*q_q)) {
        *out_q = std::numeric_limits<RealT>::quiet_NaN();
      } else {
        using MajiqInclude::quantile;
        const npy_double compute_q =
            *q_q < 0 ? npy_double{0} :
            *q_q > 1 ? npy_double{1} : npy_double{*q_q};
        *out_q = quantile(x_notnan.begin(), x_notnan_end, compute_q);
      }
    }
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
  NPY_FLOAT, NPY_DOUBLE, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace NanQuantile
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_NANQUANTILE_HPP
