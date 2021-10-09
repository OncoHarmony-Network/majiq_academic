/**
 * TTestSample.hpp
 *
 * Compute quantiles of ttest on beta distribution mixture samples
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Simultaneous sampling, testing, and summarization of mixture distributions
 * from two groups of experiments. This enables increased memory efficiency
 * (and likely small speedup as well).
 * Given J junctions, N samples, M bootstrap replicates and P psisamples
 * desired, we have:
 * a[J,N,M], b[J,N,M] --> x[J,P,N] samples from mixture distribution
 * x[J,P,N], labels[N] --> p[J,P] pvalue samples
 * p[J,P], q --> p_quantile[J] pvalue quantile
 *
 * When P > M, intermediate x[J,P,N] is the largest array in memory.
 * Simultaneous sampling and testing allows us to go straight to p[J,P] which
 * would only be largest if P > N * M.
 * This implementation goes straight to p_quantile[J], which will never be the
 * largest object.
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_TTESTSAMPLE_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_TTESTSAMPLE_HPP

#include <numpy/ndarraytypes.h>

#include <limits>
#include <vector>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/RNGPool.hpp>
#include <majiqinclude/quantile.hpp>
#include <majiqinclude/BetaMixture.hpp>
#include <majiqinclude/stats/TTest.hpp>

#include <boost/random/mersenne_twister.hpp>


namespace MajiqGufuncs {
namespace BetaMixture {
namespace TTestSample {

static char name[] = "ttest_sample";
constexpr int nin = 5;
constexpr int nout = 1;
static char signature[] = "(n,m),(n,m),(n),(q?),()->(q?)";
static char doc[] = R"pbdoc(
Obtain samples for uniform mixture of beta distributions

Parameters
----------
a, b: array[float]
    The parameters of the mixture of the experiments' beta distributions.
    core axes (n, m) correspond to (n) experiments, (m) mixture parameters.
    NaN values indicate missing observations.
    Testing performed on samples from these distributions.
labels: array[bool]
    Assignment into one of two groups for testing.
q: array[float]
    Quantiles of sampled test statistics to return
psisamples: int
    Number of sampled test statistics to take quantiles from
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
  const npy_intp dim_exp = dimensions[0];
  const npy_intp dim_mix = dimensions[1];
  const npy_intp dim_q = dimensions[2];
  // inner strides
  const npy_intp str_a_exp = steps[0];
  const npy_intp str_a_mix = steps[1];
  const npy_intp str_b_exp = steps[2];
  const npy_intp str_b_mix = steps[3];
  const npy_intp str_labels_exp = steps[4];
  const npy_intp str_q_q = steps[5];
  const npy_intp str_out_q = steps[6];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto a = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto b = CoreIt<RealT>::begin(args[1], outer_stride[1]);
  auto labels = CoreIt<bool>::begin(args[2], outer_stride[2]);
  auto q = CoreIt<RealT>::begin(args[3], outer_stride[3]);
  auto psisamples = CoreIt<int64_t>::begin(args[4], outer_stride[4]);
  auto out = CoreIt<RealT>::begin(args[5], outer_stride[5]);

  if (dim_exp < 1 || dim_mix < 1) {
    // no samples
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++out) {
      out
        .with_stride(str_out_q)
        .fill(dim_q, std::numeric_limits<RealT>::quiet_NaN());
    }
  }
  // otherwise
  // acquire random number generator
  using MajiqInclude::RNGPool;
  using boost::random::mt19937;
  RNGPool<mt19937>& rng_pool = *static_cast<RNGPool<mt19937>*>(data);
  auto rng_ptr = rng_pool.acquire();
  mt19937& gen = *rng_ptr;
  // buffer for dim_exp samples from distributions
  std::vector<RealT> x(dim_exp);
  // outer loop for sampling
  for (npy_intp i = 0; i < dim_broadcast;
      ++i, ++a, ++b, ++labels, ++q, ++psisamples, ++out) {
    // store repeated samples of p-values from test
    std::vector<RealT> pvalue_samples(*psisamples);
    for (auto& pval : pvalue_samples) {
      // sample from each distribution, fill x
      {
        auto a_exp = a.with_stride(str_a_exp);
        auto b_exp = b.with_stride(str_b_exp);
        for (npy_intp j = 0; j < dim_exp; ++j, ++a_exp, ++b_exp) {
          auto a_mix = a_exp.with_stride(str_a_mix);
          auto b_mix = b_exp.with_stride(str_b_mix);
          using MajiqInclude::BetaMixture::IsInvalid;
          using MajiqInclude::BetaMixture::_SampleMixture_unchecked;
          x[j] = IsInvalid(a_mix, b_mix, dim_mix)
            ? std::numeric_limits<RealT>::quiet_NaN()
            : _SampleMixture_unchecked(gen, a_mix, b_mix, dim_mix);
        }  // sample from each mixture distribution
      }
      // perform test on x, labels, set to pval
      using MajiqInclude::TTest::Test;
      pval = Test(x.begin(), labels.with_stride(str_labels_exp), dim_exp);
    }  // compute psisamples pvalues
    // obtain desired quantiles from pvalue_samples, set to out
    {
      auto q_q = q.with_stride(str_q_q);
      auto out_q = out.with_stride(str_out_q);
      for (npy_intp k = 0; k < dim_q; ++k, ++q_q, ++out_q) {
        using MajiqInclude::quantile;
        *out_q = quantile(pvalue_samples.begin(), pvalue_samples.end(), *q_q);
      }  // loop over quantiles
    }  // done taking quantiles of p-values
  }  // done loop over broadcast dimensions
  return;
}

constexpr int ntypes = 2;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>),
};
static char types[ntypes * (nin + nout)] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_FLOAT, NPY_BOOL, NPY_FLOAT, NPY_INT64, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_BOOL, NPY_DOUBLE, NPY_INT64, NPY_DOUBLE,
};

}  // namespace TTestSample
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_TTESTSAMPLE_HPP
