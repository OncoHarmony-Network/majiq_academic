/**
 * MeansMedian.hpp
 *
 * Calculate median of means of components of beta mixture
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_MEANSMEDIAN_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_MEANSMEDIAN_HPP

#include <numpy/ndarraytypes.h>

#include <limits>
#include <vector>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/BetaMixture.hpp>

namespace MajiqGufuncs {
namespace BetaMixture {
namespace MeansMedian {

static char name[] = "means_median";
constexpr int nin = 2;
constexpr int nout = 1;
static char signature[] = "(m),(m)->()";
static char doc[] = R"pbdoc(
Compute median of means of the beta distribution mixture

Parameters
----------
a, b: array[float]
  The parameters of the mixture of beta distributions (mixture components
  enumerated on core axis)
out: array[float]
  Output array with median(means)
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

  // pointers/iterators to data
  using MajiqGufuncs::detail::CoreIt;
  auto a = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto b = CoreIt<RealT>::begin(args[1], outer_stride[1]);
  auto out = CoreIt<RealT>::begin(args[2], outer_stride[2]);

  if (dim_mixture < 1) {
    out.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    return;
  }
  // otherwise, loop on broadcasted variables
  // make buffer to hold means that can be used for each entry
  std::vector<RealT> means(dim_mixture);
  auto means_it = means.begin();
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++a, ++b, ++out) {
    using MajiqInclude::BetaMixture::_MedianOfMeans;
    *out = _MedianOfMeans(
        a.with_stride(inner_stride_a), b.with_stride(inner_stride_b),
        means_it, dim_mixture);
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

}  // namespace MeansMedian
}  // namespace BetaMixture
}  // namespace MajiqGufuncs
#endif  // MAJIQGUFUNCS_BETAMIXTURE_MEANSMEDIAN_HPP
