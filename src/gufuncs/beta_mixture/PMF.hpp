/**
 * PMF.hpp
 *
 * Implementation of outer loop for beta mixture PMF ufunc
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_PMF_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_PMF_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/BetaMixture.hpp>


namespace MajiqGufuncs {
namespace BetaMixture {
namespace PMF {

static char name[] = "pmf";
constexpr int nin = 3;
constexpr int nout = 1;
static char signature[] = "(m),(m),(b)->(b)";
static char doc[] = R"pbdoc(
Compute PMF of mixture of beta distributions

Parameters
----------
a, b: array[float]
    The parameters of the mixture of beta distributions (mixture components
    enumerated on core axis)
bins_dummy: array[float]
    sentinel array used for determining number of desired bins
)pbdoc";

template <typename RealT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides for broadcasting
  const npy_intp* stride_broadcast = steps;
  steps += nin + nout;

  // NOTE: we skip steps[2], args[2], stride_broadcast[2] because those
  // correspond to dummy variable not used in computation besides letting numpy
  // know the desired size of out

  // core dimensions
  const npy_intp dim_mixture = dimensions[0];
  const npy_intp dim_bins = dimensions[1];
  // strides on core dimensions
  const npy_intp inner_stride_a = steps[0];
  const npy_intp inner_stride_b = steps[1];
  const npy_intp inner_stride_out = steps[3];

  // pointers/iterators to data
  using MajiqGufuncs::detail::CoreIt;
  auto a = CoreIt<RealT>::begin(args[0], stride_broadcast[0]);
  auto b = CoreIt<RealT>::begin(args[1], stride_broadcast[1]);
  auto out = CoreIt<RealT>::begin(args[3], stride_broadcast[3]);

  if (dim_mixture < 1) {
    // for each iteration of broadcasting, fill core dimension of out with nan
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++out) {
      out.with_stride(inner_stride_out)
        .fill(dim_bins, std::numeric_limits<RealT>::quiet_NaN());
    }
    return;
  }
  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++a, ++b, ++out) {
    using MajiqInclude::BetaMixture::_PMF;
    _PMF(a.with_stride(inner_stride_a),
        b.with_stride(inner_stride_b),
        out.with_stride(inner_stride_out),
        dim_mixture, dim_bins);
  }
  return;
}

constexpr int ntypes = 3;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>),
};
static char types[ntypes * (nin + nout)] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
  NPY_FLOAT, NPY_FLOAT, NPY_DOUBLE, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace PMF
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif
