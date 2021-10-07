/**
 * Approximation.hpp
 *
 * Calculate first approximation of beta mixture
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_APPROXIMATION_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_APPROXIMATION_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include <gufuncs/helpers.hpp>
#include "BetaMixture.hpp"

namespace MajiqGufuncs {
namespace BetaMixture {
namespace Approximation {

static char name[] = "approximation";
constexpr int nin = 2;
constexpr int nout = 1;
static char signature[] = "(m),(m)->(2)";
static char doc[] = R"pbdoc(
Compute beta distribution parameters matching beta mixture mean, variance

Parameters
----------
a, b: array[float]
  The parameters of the mixture of beta distributions (mixture components
  enumerated on core axis)
out: array[float]
  Output array with a0 = out[..., 0], b0 = out[..., 1] parameters of beta
  distribution matching mean, variance of mixture
)pbdoc";

template <typename RealT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  const npy_intp dim_broadcast = *dimensions++;
  const npy_intp* outer_stride = steps;
  steps += nin + nout;

  // core dimensions
  const npy_intp dim_mixture = dimensions[0];
  // inner strides
  const npy_intp inner_stride_a = steps[0];
  const npy_intp inner_stride_b = steps[1];
  const npy_intp inner_stride_out = steps[2];

  // pointers/iterators to data
  auto a = detail::CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto b = detail::CoreIt<RealT>::begin(args[1], outer_stride[1]);
  auto a0 = detail::CoreIt<RealT>::begin(args[2], outer_stride[2]);
  auto b0 = (
      1 + a0.with_stride(inner_stride_out)).with_stride(outer_stride[2]);

  if (dim_mixture < 1) {
    a0.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    b0.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
  } else if (dim_mixture == 1) {
    // identity -- copy values
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++a, ++b, ++a0, ++b0) {
      *a0 = *a;
      *b0 = *b;
    }
  } else {
    // otherwise, outer loop on broadcasted variables
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++a, ++b, ++a0, ++b0) {
      auto approximation = _BetaApproximation(a.with_stride(inner_stride_a),
          b.with_stride(inner_stride_b), dim_mixture);
      *a0 = approximation.first;
      *b0 = approximation.second;
    }
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
  NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace Approximation
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_APPROXIMATION_HPP
