/**
 * Sample.hpp
 *
 * Sample from beta distribution
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_SAMPLING_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_SAMPLING_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include "GlobalRNGPool.hpp"

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/BetaMixture.hpp>



namespace MajiqGufuncs {
namespace BetaMixture {
namespace Sample {

static char name[] = "sample";
constexpr int nin = 2;
constexpr int nout = 1;
static char signature[] = "(m),(m)->()";
static char doc[] = R"pbdoc(
Obtain samples for uniform mixture of beta distributions

Parameters
----------
a, b: array[float]
    The parameters of the mixture of beta distributions (mixture components
    enumerated on core axis)
out: array[float]
    The output array with corect size that will be filled in

Notes
-----
Random number generation can be controlled by rng_seed() and rng_resize().
Random number generation in this module is threadsafe; rng_resize() sets the
number of random number generators that can be made available to the different
random number generators.
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
  const npy_intp dim_mixture = dimensions[0];
  // inner strides
  const npy_intp inner_stride_a = steps[0];
  const npy_intp inner_stride_b = steps[1];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto a = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto b = CoreIt<RealT>::begin(args[1], outer_stride[1]);
  auto out = CoreIt<RealT>::begin(args[2], outer_stride[2]);

  if (dim_mixture < 1) {
    // no samples to process, so must be nan
    out.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    return;
  }
  // otherwise

  // acquire ownership random number generator
  // NOTE: need to keep pointer in scope to maintain ownership
  // (i.e. don't replace with *global_rng_pool.acquire())
  auto gen_ptr = global_rng_pool.acquire();
  auto& gen = *gen_ptr;

  // outer loop for sampling
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++a, ++b, ++out) {
    auto inner_a = a.with_stride(inner_stride_a);
    auto inner_b = b.with_stride(inner_stride_b);
    using MajiqInclude::BetaMixture::IsInvalid;
    using MajiqInclude::BetaMixture::_SampleMixture_unchecked;
    *out = IsInvalid(inner_a, inner_b, dim_mixture)
      ? std::numeric_limits<RealT>::quiet_NaN()
      : _SampleMixture_unchecked(gen, inner_a, inner_b, dim_mixture);
  }
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

}  // namespace Sample
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_SAMPLING_HPP