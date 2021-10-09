/**
 * Moments.hpp
 *
 * Calculate first moments of beta mixture
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_MOMENTS_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_MOMENTS_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/BetaMixture.hpp>

namespace MajiqGufuncs {
namespace BetaMixture {
namespace Moments {

static char name[] = "moments";
constexpr int nin = 2;
constexpr int nout = 1;
static char signature[] = "(m),(m)->(2)";
static char doc[] = R"pbdoc(
Compute mean and variance of beta distribution mixture

Parameters
----------
a, b: array[float]
  The parameters of the mixture of beta distributions (mixture components
  enumerated on core axis)
out: array[float]
  Output array with mean = out[..., 0], var = out[..., 1]
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
  using MajiqGufuncs::detail::CoreIt;
  auto a = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto b = CoreIt<RealT>::begin(args[1], outer_stride[1]);
  auto mean = CoreIt<RealT>::begin(args[2], outer_stride[2]);
  auto var = (1 + mean.with_stride(inner_stride_out))
    .with_stride(outer_stride[2]);

  if (dim_mixture < 1) {
    mean.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    var.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    return;
  }
  // otherwise, outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++a, ++b, ++mean, ++var) {
    using MajiqInclude::BetaMixture::_Moments;
    auto moments = _Moments(a.with_stride(inner_stride_a),
        b.with_stride(inner_stride_b), dim_mixture);
    *mean = moments.mean;
    *var = moments.variance;
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

}  // namespace Moments
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_MOMENTS_HPP
