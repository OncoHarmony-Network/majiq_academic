/**
 * PDF.hpp
 *
 * Implementation of outer loop for beta mixture PDF ufunc
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_PDF_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_PDF_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include <gufuncs/helpers.hpp>
#include "BetaMixture.hpp"


namespace MajiqGufuncs {
namespace BetaMixture {
namespace PDF {

static char name[] = "pdf";
constexpr int nin = 3;
constexpr int nout = 1;
static char signature[] = "(),(m),(m)->()";
static char doc[] = R"pbdoc(
Compute PDF of mixture of beta distributions

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
  const npy_intp stride_x = *steps++;
  const npy_intp stride_a = *steps++;
  const npy_intp stride_b = *steps++;
  const npy_intp stride_out = *steps++;
  // core dimensions
  const npy_intp dim_mixture = dimensions[0];
  // inner strides
  const npy_intp inner_stride_a = steps[0];
  const npy_intp inner_stride_b = steps[1];
  // pointers/iterators to data
  auto x = detail::CoreIt<RealT>::begin(args[0], stride_x);
  auto a = detail::CoreIt<RealT>::begin(args[1], inner_stride_a);
  auto b = detail::CoreIt<RealT>::begin(args[2], inner_stride_b);
  auto out = detail::CoreIt<RealT>::begin(args[3], stride_out);

  if (dim_mixture < 1) {
    // if there are no distributions, PDF is NaN everywhere
    if (dim_broadcast > 1) {
      std::fill(
          out, out + dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    } else if (dim_broadcast == 1) {
      out[0] = std::numeric_limits<RealT>::quiet_NaN();
    }
    return;
  }
  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++out,
      a.apply_stride(stride_a), b.apply_stride(stride_b)) {
    *out = _PDF(*x, a, b, dim_mixture);
  }
  return;
}

constexpr int ntypes = 2;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>),
};
static char types[
  ntypes * (nin + nout)
] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace PDF
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif
