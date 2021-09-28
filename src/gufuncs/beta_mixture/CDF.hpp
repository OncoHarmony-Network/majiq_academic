/**
 * CDF.hpp
 *
 * Implementation of outer loop for beta mixture CDF ufunc
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_CDF_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_CDF_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include <gufuncs/helpers.hpp>
#include "BetaMixture.hpp"


namespace MajiqGufuncs {
namespace BetaMixture {
namespace CDF {

static char name[] = "cdf";
constexpr int nin = 3;
constexpr int nout = 1;
static char signature[] = "(),(m),(m)->()";
static char doc[] = R"pbdoc(
Compute CDF of mixture of beta distributions

Parameters
----------
x: array[float]
    Points at which to evaluate cdf
a, b: array[float]
    The parameters of the mixture of beta distributions (mixture components
    enumerated on core axis)
out: array[float]
    The output array with correct size that will be filled in
)pbdoc";

template <typename RealT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides on each variable for outer loop
  const npy_intp* stride_broadcast = steps;
  steps += nin + nout;

  // core dimensions
  const npy_intp dim_mixture = dimensions[0];
  // inner strides
  const npy_intp inner_stride_a = steps[0];
  const npy_intp inner_stride_b = steps[1];

  // pointers/iterators to data
  auto x = detail::CoreIt<RealT>::begin(args[0], stride_broadcast[0]);
  auto a = detail::CoreIt<RealT>::begin(args[1], stride_broadcast[1]);
  auto b = detail::CoreIt<RealT>::begin(args[2], stride_broadcast[2]);
  auto out = detail::CoreIt<RealT>::begin(args[3], stride_broadcast[3]);

  if (dim_mixture < 1) {
    out.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    return;
  }
  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++a, ++b, ++out) {
    *out = _CDF(*x, a.with_stride(inner_stride_a),
        b.with_stride(inner_stride_b), dim_mixture);
  }
  return;
}

constexpr int ntypes = 2;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>),
};
static char types[ntypes * (nin + nout)] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace CDF
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif
