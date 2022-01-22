/**
 * StochasticPMF.hpp
 *
 * Approximation of approximate PMF
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_STOCHASTICPMF_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_STOCHASTICPMF_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include "GlobalRNGPool.hpp"

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/BetaMixture.hpp>

namespace MajiqGufuncs {
namespace BetaMixture {
namespace StochasticPMF {

static char name[] = "stochastic_pmf";
constexpr int nin = 4;
constexpr int nout = 1;
static char signature[] = "(m),(m),(b),()->(b)";
static char doc[] = R"pbdoc(
stochastic_pmf(a, b, bins_dummy, psisamples)

Average histogram of psisamples draws from each beta distribution in a mixture

Notes
-----
Will ignore invalid distributions (i.e. total number of samples will equal
psisamples * number_valid_distributions).
This is different than other beta mixture functions.

Parameters
----------
a, b: array[float]
    The parameters of the mixture of beta distributions (mixture components
    enumerated on core axis)
bins_dummy: array[float]
    sentinel array used for determining number of desired bins
psisamples: int
    Number of test statistics from distribution samples to take quantiles from
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
  auto psisamples = CoreIt<int64_t>::begin(args[3], stride_broadcast[3]);
  auto out = CoreIt<int64_t>::begin(args[4], stride_broadcast[4]);

  if (dim_bins < 1) { return; }
  if (dim_mixture < 1) {
    // for each iteration of broadcasting fill output with 0
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++out) {
      // fill output with 0
      out.with_stride(inner_stride_out).fill(dim_bins, int64_t{0});
    }
    return;
  }

  // acquire ownership random number generator
  // NOTE: need to keep pointer in scope to maintain ownership
  // (i.e. don't replace with *global_rng_pool.acquire())
  RNGPtrT gen_ptr = global_rng_pool.acquire();

  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++a, ++b, ++psisamples, ++out) {
    // get vectors over core dimensions
    auto a_mix = a.with_stride(inner_stride_a);
    auto b_mix = b.with_stride(inner_stride_b);
    auto out_bins = out.with_stride(inner_stride_out);
    out_bins.fill(dim_bins, int64_t{0});
    // how many psisamples?
    if (*psisamples < 1) { continue; }  // no samples if no psisamples
    // loop over each mixture component
    for (npy_intp j = 0; j < dim_mixture; ++j, ++a_mix, ++b_mix) {
      // skip component if invalid
      using MajiqInclude::BetaMixture::IsInvalidComponent;
      if (IsInvalidComponent(*a_mix, *b_mix)) { continue; }
      // otherwise take psisamples from this component
      using MajiqInclude::BetaMixture::_Sample_unchecked;
      for (int64_t k = 0; k < *psisamples; ++k) {
        // sample from mixture distribution
        const auto x = _Sample_unchecked(*gen_ptr, *a_mix, *b_mix);
        // bin which it belongs to
        auto idx = static_cast<npy_intp>(std::floor(dim_bins * x));
        if (idx < 0) {
          idx = 0;
        } else if (idx >= dim_mixture) {
          idx = dim_mixture - 1;
        }
        ++out_bins[idx];
      }
    }
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
  NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_INT64, NPY_INT64,
  NPY_FLOAT, NPY_FLOAT, NPY_DOUBLE, NPY_INT64, NPY_INT64,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_INT64, NPY_INT64,
};

}  // namespace StochasticPMF
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_STOCHASTICPMF_HPP
