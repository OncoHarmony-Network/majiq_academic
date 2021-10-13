/**
 * StatsSample.hpp
 *
 * Compute quantiles of multiple test statistics
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

#ifndef MAJIQGUFUNCS_BETAMIXTURE_STATSSAMPLE_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_STATSSAMPLE_HPP

#include <numpy/ndarraytypes.h>

#include <algorithm>
#include <array>
#include <limits>
#include <set>
#include <vector>

#include "GlobalRNGPool.hpp"

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/quantile.hpp>
#include <majiqinclude/BetaMixture.hpp>
#include <majiqinclude/stats/InfoScore.hpp>
#include <majiqinclude/stats/MannWhitney.hpp>
#include <majiqinclude/stats/TNOM.hpp>
#include <majiqinclude/stats/TTest.hpp>


namespace MajiqGufuncs {
namespace BetaMixture {
namespace StatsSample {

enum class HetStats : int64_t {
  TTest,
  MannWhitney,
  TNOM,
  InfoScore,
  Count  // this one has value equal to number of valid states
};

static char name[] = "stats_sample";
constexpr int nin = 6;
constexpr int nout = 1;
static char signature[] = "(n,m),(n,m),(n),(q?),(),(stat?)->(stat?,q?)";
static char doc[] = R"pbdoc(
Obtain samples for uniform mixture of beta distributions

Parameters
----------
a, b: array[float] (n, m)
    The parameters of the mixture of the experiments' beta distributions.
    core axes (n, m) correspond to (n) experiments, (m) mixture parameters.
    NaN values indicate missing observations.
    Testing performed on samples from these distributions.
labels: array[bool] (n)
    Assignment into one of two groups for testing.
q: array[float] (q?)
    Quantiles of sampled test statistics to return
psisamples: int ()
    Number of sampled test statistics to take quantiles from
stats: array[int] (stat?)
    Statistics to use. 0 = ttest, 1 = MannWhitney, 2 = TNOM, 3 = InfoScore
out: array[float] (stat?,q?)
    The output array with corect size that will be filled in
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
  const npy_intp dim_stat = dimensions[3];
  // inner strides
  const npy_intp str_a_exp = steps[0];
  const npy_intp str_a_mix = steps[1];
  const npy_intp str_b_exp = steps[2];
  const npy_intp str_b_mix = steps[3];
  const npy_intp str_labels_exp = steps[4];
  const npy_intp str_q_q = steps[5];
  const npy_intp str_stats_stat = steps[6];
  const npy_intp str_out_stat = steps[7];
  const npy_intp str_out_q = steps[8];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto a = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto b = CoreIt<RealT>::begin(args[1], outer_stride[1]);
  auto labels = CoreIt<bool>::begin(args[2], outer_stride[2]);
  auto q = CoreIt<RealT>::begin(args[3], outer_stride[3]);
  auto psisamples = CoreIt<int64_t>::begin(args[4], outer_stride[4]);
  auto stats = CoreIt<int64_t>::begin(args[5], outer_stride[5]);
  auto out = CoreIt<double>::begin(args[6], outer_stride[6]);

  // if any of the core dimensions are empty, output will be trivial
  if (dim_q < 1 || dim_stat < 1) {
    // output will be empty!
    return;
  }
  if (dim_exp < 1 || dim_mix < 1) {
    // no samples
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++out) {
      auto out_stat = out.with_stride(str_out_stat);
      for (npy_intp z = 0; z < dim_stat; ++z, ++out_stat) {
        out_stat
          .with_stride(str_out_q)
          .fill(dim_q, std::numeric_limits<RealT>::quiet_NaN());
      }
    }
  }
  // otherwise

  // acquire ownership random number generator
  // NOTE: need to keep pointer in scope to maintain ownership
  // (i.e. don't replace with *global_rng_pool.acquire())
  auto gen_ptr = global_rng_pool.acquire();
  auto& gen = *gen_ptr;

  // caching test results for non-ttest statistics
  MajiqInclude::MannWhitney::MannWhitneyCache mannwhitney;
  MajiqInclude::TNOM::TNOMCache tnom;
  MajiqInclude::InfoScore::InfoScoreCache infoscore;

  // buffer for dim_exp samples from respective distributions
  std::vector<RealT> x(dim_exp);

  // buffer to indicate if distribution is invalid for a given set of
  // parameters (check only for first sample)
  std::vector<bool> is_invalid(dim_exp);

  // buffer to indicate how to sort x
  std::vector<size_t> sortx(dim_exp);
  std::iota(sortx.begin(), sortx.end(), size_t{0});

  // buffer for sampled statistics
  std::array<std::vector<double>, static_cast<size_t>(HetStats::Count)> pvalue_samples;

  // buffer for quantiles of statistics from samples
  std::array<std::vector<double>, static_cast<size_t>(HetStats::Count)> pvalue_quantiles;
  std::for_each(pvalue_quantiles.begin(), pvalue_quantiles.end(),
      [dim_q](std::vector<double>& q) { q.resize(dim_q); });

  // outer loop
  for (npy_intp i = 0; i < dim_broadcast;
      ++i, ++a, ++b, ++labels, ++q, ++psisamples, ++stats, ++out) {
    auto stats_stat = stats.with_stride(str_stats_stat);
    const int64_t n_psisamples{*psisamples};

    // identify which statistics to compute.
    std::set<HetStats> compute_stat{};
    if (n_psisamples > 0) {
      // only compute if we plan to do any psisamples
      for (npy_intp z = 0; z < dim_stat; ++z) {
        const auto& stat_z = stats_stat[z];
        if (stat_z >= 0 && stat_z < static_cast<int64_t>(HetStats::Count)) {
          compute_stat.insert(static_cast<HetStats>(stat_z));
        }
      }
    }

    // if there are any statistics to compute, do so
    if (!compute_stat.empty()) {
      // prepare output buffers for requested statistics per psisample
      for (const auto& stat : compute_stat) {
        pvalue_samples[static_cast<int64_t>(stat)].resize(n_psisamples);
      }
      // for each psisample, sample x, update psisamples
      // we check if compute_stat is empty because we remove if returns nan
      for (int64_t s = 0; s < n_psisamples && !compute_stat.empty(); ++s) {
        // sample x
        auto a_exp = a.with_stride(str_a_exp);
        auto b_exp = b.with_stride(str_b_exp);
        for (npy_intp j = 0; j < dim_exp; ++j, ++a_exp, ++b_exp) {
          auto a_mix = a_exp.with_stride(str_a_mix);
          auto b_mix = b_exp.with_stride(str_b_mix);
          if (s == 0) {
            // first psisample, we check whether or not it's valid
            using MajiqInclude::BetaMixture::IsInvalid;
            is_invalid[j] = IsInvalid(a_mix, b_mix, dim_mix);
          }
          using MajiqInclude::BetaMixture::_SampleMixture_unchecked;
          x[j] = is_invalid[j] ? std::numeric_limits<RealT>::quiet_NaN()
            : _SampleMixture_unchecked(gen, a_mix, b_mix, dim_mix);
        }  // done sampling x
        // sort x if we want a non-ttest statistic
        if (std::find_if(compute_stat.begin(), compute_stat.end(),
              [](const HetStats& u) -> bool { return u != HetStats::TTest; })
            != compute_stat.end()) {
          // we want a statistic that isn't t-test. This will require sorting.
          // numpy sorts nans to end and that's what we want. To do this with
          // C++ STL, we partition on nan vs not, then sort the non-nan portion.
          auto first_nan = std::partition(sortx.begin(), sortx.end(),
              [&x](size_t i) { return !std::isnan(x[i]); });
          std::sort(sortx.begin(), first_nan,
              [&x](size_t i, size_t j) { return x[i] < x[j]; });
        }
        // perform desired tests
        std::set<HetStats> delete_stat;
        for (const auto& stat : compute_stat) {
          const auto stat_idx = static_cast<int64_t>(stat);
          switch (stat) {
            case HetStats::TTest:
              pvalue_samples[stat_idx][s] = MajiqInclude::TTest::Test(
                  x.begin(), labels.with_stride(str_labels_exp), dim_exp);
              if (std::isnan(pvalue_samples[stat_idx][s])) {
                delete_stat.insert(stat);
              }
              break;
            case HetStats::MannWhitney:
              pvalue_samples[stat_idx][s] = MajiqInclude::MannWhitney::Test(
                  mannwhitney, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              if (std::isnan(pvalue_samples[stat_idx][s])) {
                delete_stat.insert(stat);
              }
              break;
            case HetStats::TNOM:
              pvalue_samples[stat_idx][s] = MajiqInclude::TNOM::Test(
                  tnom, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              if (std::isnan(pvalue_samples[stat_idx][s])) {
                delete_stat.insert(stat);
              }
              break;
            case HetStats::InfoScore:
              pvalue_samples[stat_idx][s] = MajiqInclude::InfoScore::Test(
                  infoscore, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              if (std::isnan(pvalue_samples[stat_idx][s])) {
                delete_stat.insert(stat);
              }
              break;
            default:
              break;
          }  // switch over possible tests
        }  // for loop over stats
        for (auto stat : delete_stat) {
          // don't compute stats that have returned nan on previous psisamples
          compute_stat.erase(stat);
        }
      }  // done doing psisamples (sampling and testing per psisample)

      // compute quantiles for stats
      for (auto stat : compute_stat) {
        auto q_q = q.with_stride(str_q_q);
        auto& in_pvalues = pvalue_samples[static_cast<int64_t>(stat)];
        auto& out_quantiles = pvalue_quantiles[static_cast<int64_t>(stat)];
        for (npy_intp k = 0; k < dim_q; ++k, ++q_q) {
          using MajiqInclude::quantile;
          out_quantiles[k] = quantile(
              in_pvalues.begin(), in_pvalues.end(), *q_q);
        }  // loop over quantiles for given stat
      }  // loop over stats we are computing
    }  // done if we were computing stats

    // write pvalue quantiles to output
    auto out_stat = out.with_stride(str_out_stat);
    for (npy_intp z = 0; z < dim_stat; ++z, ++out_stat) {
      const auto& stat_z = stats_stat[z];
      if (compute_stat.count(static_cast<HetStats>(stat_z)) > 0) {
        std::copy(
            pvalue_quantiles[stat_z].begin(), pvalue_quantiles[stat_z].end(),
            out_stat.with_stride(str_out_q));
      } else {
        // wasn't computed, so NaN
        out_stat
          .with_stride(str_out_q)
          .fill(dim_q, std::numeric_limits<RealT>::quiet_NaN());
      }  // copy over computed quantiles or nan
    }  // loop over requested statistics
  }  // outer loop
  return;
}

constexpr int ntypes = 2;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>),
};
static char types[ntypes * (nin + nout)] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_FLOAT, NPY_BOOL, NPY_FLOAT, NPY_INT64, NPY_INT64, NPY_DOUBLE,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_BOOL, NPY_DOUBLE, NPY_INT64, NPY_INT64, NPY_DOUBLE,
};

}  // namespace StatsSample
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_STATSSAMPLE_HPP
