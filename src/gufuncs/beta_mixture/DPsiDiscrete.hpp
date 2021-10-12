/**
 * DPsiDiscrete.hpp
 *
 * Given discretized log prior and beta mixture for psi1, psi2, compute log
 * posterior
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_DPSIDISCRETE_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_DPSIDISCRETE_HPP

#include <numpy/ndarraytypes.h>

#include <limits>
#include <vector>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/BetaMixture.hpp>
#include <majiqinclude/DPsiDiscrete.hpp>


namespace MajiqGufuncs {
namespace BetaMixture {
namespace DPsiDiscrete {

static char name[] = "dpsi_discrete";
constexpr int nin = 5;
constexpr int nout = 1;
static char signature[] = "(m1),(m1),(m2),(m2),(b)->(b)";
static char doc[] = R"pbdoc(
Compute discrete log posterior of difference between distributions given prior

Parameters
----------
a1, b1, a2, b2: array[float]
    The parameters of the mixtures of beta distributions (mixture components
    enumerated on core axis)
logprior: array[float]
    log prior for deltapsi distribution. The number of bins must be an even
    number (otherwise will return all nan).
)pbdoc";

template <typename RealT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides for broadcasting
  const npy_intp* stride_broadcast = steps;
  steps += nin + nout;

  // core dimensions
  const npy_intp dim_mixture1 = dimensions[0];
  const npy_intp dim_mixture2 = dimensions[1];
  const npy_intp dim_bins = dimensions[2];
  // strides on core dimensions
  const npy_intp inner_stride_a1 = steps[0];
  const npy_intp inner_stride_b1 = steps[1];
  const npy_intp inner_stride_a2 = steps[2];
  const npy_intp inner_stride_b2 = steps[3];
  const npy_intp inner_stride_prior = steps[4];
  const npy_intp inner_stride_out = steps[5];

  // pointers/iterators to data
  using MajiqGufuncs::detail::CoreIt;
  auto a1 = CoreIt<RealT>::begin(args[0], stride_broadcast[0]);
  auto b1 = CoreIt<RealT>::begin(args[1], stride_broadcast[1]);
  auto a2 = CoreIt<RealT>::begin(args[2], stride_broadcast[2]);
  auto b2 = CoreIt<RealT>::begin(args[3], stride_broadcast[3]);
  auto logprior = CoreIt<RealT>::begin(args[4], stride_broadcast[4]);
  auto out = CoreIt<RealT>::begin(args[5], stride_broadcast[5]);

  if (dim_mixture1 < 1 || dim_mixture2 < 1 || dim_bins % 2 != 0) {
    // for each iteration of broadcasting, fill core dimension of out with nan
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++out) {
      out.with_stride(inner_stride_out)
        .fill(dim_bins, std::numeric_limits<RealT>::quiet_NaN());
    }
    return;
  }
  // outer loop on broadcasted variables
  int64_t n_psibins = dim_bins / 2;  // bins for psi distribution
  // buffers for psi1, psi2
  std::vector<RealT> logp_psi1(n_psibins);
  std::vector<RealT> logp_psi2(n_psibins);
  for (npy_intp i = 0; i < dim_broadcast;
      ++i, ++a1, ++b1, ++a2, ++b2, ++logprior, ++out) {
    // get buffer with stride for inner loop
    auto a1_mix = a1.with_stride(inner_stride_a1);
    auto a2_mix = a2.with_stride(inner_stride_a2);
    auto b1_mix = b1.with_stride(inner_stride_b1);
    auto b2_mix = b2.with_stride(inner_stride_b2);
    auto out_bins = out.with_stride(inner_stride_out);
    // if invalid beta mixtures for psi1 or psi2, make nan
    using MajiqInclude::BetaMixture::IsInvalid;
    if (IsInvalid(a1_mix, b1_mix, dim_mixture1)
        || IsInvalid(a2_mix, b2_mix, dim_mixture2)) {
      out_bins.fill(dim_bins, std::numeric_limits<RealT>::quiet_NaN());
    } else {
      // compute PMFs
      using MajiqInclude::BetaMixture::_LogPMF;
      // compute log posterior for psi1, psi2
      _LogPMF(a1.with_stride(inner_stride_a1), b1.with_stride(inner_stride_b1),
          logp_psi1.begin(), dim_mixture1, n_psibins);
      _LogPMF(a2.with_stride(inner_stride_a2), b2.with_stride(inner_stride_b2),
          logp_psi2.begin(), dim_mixture2, n_psibins);
      // get unnormalized log distribution convolving/interpolating for dpsi
      // (without dpsi prior)
      using MajiqInclude::DPsiDiscrete::DPsiIndependent;
      DPsiIndependent(
          logp_psi1.begin(), logp_psi2.begin(), n_psibins, out_bins);
      // apply prior, renormalize
      using MajiqInclude::DPsiDiscrete::DPsiPosterior;
      DPsiPosterior(
          out_bins, logprior.with_stride(inner_stride_prior), dim_bins);
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
  NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace DPsiDiscrete
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_DPSIDISCRETE_HPP
