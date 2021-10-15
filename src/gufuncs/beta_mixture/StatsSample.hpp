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
#include <limits>
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
constexpr int nout = 2;
static char signature[] = "(n,m),(n,m),(n),(q?),(),(stat?)->(stat?),(stat?,q?)";
static char doc[] = R"pbdoc(
stats_sample(a, b, labels, q, psisamples, stats, [out1, [out2]], /, ...)

Obtain samples for uniform mixture of beta distributions

Parameters
----------
a, b: array[float] (n, m)
    The parameters of the mixture of the experiments' beta distributions.
    core axes (n, m) correspond to (n) experiments, (m) mixture parameters.
    NaN values indicate missing observations.
    Testing performed on means and samples from these distributions.
labels: array[bool] (n)
    Assignment into one of two groups for testing.
q: array[float] (q?)
    Quantiles of sampled test statistics to return
psisamples: int ()
    Number of test statistics from distribution samples to take quantiles from
stats: array[int] (stat?)
    Statistics to use. 0 = ttest, 1 = MannWhitney, 2 = TNOM, 3 = InfoScore
out=(out1, out2): Tuple[array[float], array[float]] (stat?),(stat?,q?)
    out1: pvalues on distribution means for requested statistics
    out2: pvalue quantiles on distribution samples for requested statistics
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
  const npy_intp str_outmean_stat = steps[7];
  const npy_intp str_out_stat = steps[8];
  const npy_intp str_out_q = steps[9];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto a = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto b = CoreIt<RealT>::begin(args[1], outer_stride[1]);
  auto labels = CoreIt<bool>::begin(args[2], outer_stride[2]);
  auto q = CoreIt<RealT>::begin(args[3], outer_stride[3]);
  auto psisamples = CoreIt<int64_t>::begin(args[4], outer_stride[4]);
  auto stats = CoreIt<int64_t>::begin(args[5], outer_stride[5]);
  auto outmean = CoreIt<double>::begin(args[6], outer_stride[6]);
  auto out = CoreIt<double>::begin(args[7], outer_stride[7]);

  // if any of the core dimensions are empty, output will be trivial
  if (dim_stat < 1) {
    // output will be empty!
    return;
  }
  if (dim_exp < 1 || dim_mix < 1) {
    // no samples
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++outmean, ++out) {
      outmean
        .with_stride(str_out_stat)
        .fill(dim_stat, std::numeric_limits<double>::quiet_NaN());
      auto out_stat = out.with_stride(str_out_stat);
      for (npy_intp z = 0; z < dim_stat; ++z, ++out_stat) {
        out_stat
          .with_stride(str_out_q)
          .fill(dim_q, std::numeric_limits<double>::quiet_NaN());
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
  // function to update sortx
  auto update_sortx = [&x, &sortx]() {
    // numpy sorts nans to end and that's what we want. To do this with
    // C++ STL, we partition on nan vs not, then sort the non-nan portion.
    auto first_nan = std::partition(sortx.begin(), sortx.end(),
        [&x](size_t i) { return !std::isnan(x[i]); });
    std::sort(sortx.begin(), first_nan,
        [&x](size_t i, size_t j) { return x[i] < x[j]; });
  };

  // should we bother testing this statistic? determined by testing on means
  std::vector<bool> stat_valid(dim_stat);
  std::vector<std::vector<double>> pvalue_samples(dim_stat);
  std::vector<std::vector<double>> pvalue_quantiles(dim_stat);
  std::for_each(pvalue_quantiles.begin(), pvalue_quantiles.end(),
      [dim_q](std::vector<double>& q) { q.resize(dim_q); });

  // outer loop
  for (npy_intp i = 0; i < dim_broadcast;
      ++i, ++a, ++b, ++labels, ++q, ++psisamples, ++stats, ++outmean, ++out) {
    // identify which statistics are available
    auto stats_stat = stats.with_stride(str_stats_stat);
    npy_intp ct_valid = 0;
    for (npy_intp z = 0; z < dim_stat; ++z) {
      const auto& stat_z = stats_stat[z];
      stat_valid[z] = (stat_z >= 0)
        && (stat_z < static_cast<int64_t>(HetStats::Count));
      if (stat_valid[z]) {
        ++ct_valid;
      }
    }

    // if there are any statistics to try, do so on means
    if (ct_valid > 0) {
      auto a_exp = a.with_stride(str_a_exp);
      auto b_exp = b.with_stride(str_b_exp);
      for (npy_intp j = 0; j < dim_exp; ++j, ++a_exp, ++b_exp) {
        x[j] = MajiqInclude::BetaMixture::_Mean(
            a_exp.with_stride(str_a_mix), b_exp.with_stride(str_b_mix),
            dim_mix);
        is_invalid[j] = std::isnan(x[j]);
      }  // fill x with distribution means, mark which experiments invalid
      bool not_sorted = true;
      // perform desired tests
      auto outmean_stat = outmean.with_stride(str_outmean_stat);
      for (npy_intp z = 0; z < dim_stat; ++z, ++outmean_stat) {
        if (!stat_valid[z]) {
          *outmean_stat = std::numeric_limits<double>::quiet_NaN();
        } else {
          // use correct statistical function, sorting data once if necessary
          switch (static_cast<HetStats>(stats_stat[z])) {
            case HetStats::TTest:
              *outmean_stat = MajiqInclude::TTest::Test(
                  x.begin(), labels.with_stride(str_labels_exp), dim_exp);
              break;
            case HetStats::MannWhitney:
              if (not_sorted) { update_sortx(); not_sorted = false; }
              *outmean_stat = MajiqInclude::MannWhitney::Test(
                  mannwhitney, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              break;
            case HetStats::TNOM:
              if (not_sorted) { update_sortx(); not_sorted = false; }
              *outmean_stat = MajiqInclude::TNOM::Test(
                  tnom, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              break;
            case HetStats::InfoScore:
              if (not_sorted) { update_sortx(); not_sorted = false; }
              *outmean_stat = MajiqInclude::InfoScore::Test(
                  infoscore, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              break;
            default:
              break;
          }
          // if result was nan, we know not to do it in the future
          if (std::isnan(*outmean_stat)) {
            stat_valid[z] = false;
            --ct_valid;
          }
        }  // ifelse on if statistic was available, now also if non-nan pvalue
      }  // done loop over requested statistics
    }  // if any statistics to try, run distribution means

    // if any statistics had non-nan values, do psisamples on them
    if (ct_valid > 0) {
      const auto n_psisamples{*psisamples};
      // prepare output buffers for requested number of psisamples
      for (npy_intp z = 0; z < dim_stat; ++z) {
        if (stat_valid[z]) {
          pvalue_samples[z].resize(n_psisamples >= 0 ? n_psisamples : 0);
        }
      }
      // for each psisample, sample x, update psisamples
      for (int64_t s = 0; s < n_psisamples && ct_valid > 0; ++s) {
        // sample x
        auto a_exp = a.with_stride(str_a_exp);
        auto b_exp = b.with_stride(str_b_exp);
        for (npy_intp j = 0; j < dim_exp; ++j, ++a_exp, ++b_exp) {
          if (!is_invalid[j]) {
            auto a_mix = a_exp.with_stride(str_a_mix);
            auto b_mix = b_exp.with_stride(str_b_mix);
            using MajiqInclude::BetaMixture::_SampleMixture_unchecked;
            x[j] = _SampleMixture_unchecked(gen, a_mix, b_mix, dim_mix);
          }
        }  // done sampling x
        bool not_sorted = true;
        // statistics over samples
        for (npy_intp z = 0; z < dim_stat; ++z) {
          if (!stat_valid[z]) { continue; }
          // use correct statistical function, sorting data once if necessary
          switch (static_cast<HetStats>(stats_stat[z])) {
            case HetStats::TTest:
              pvalue_samples[z][s] = MajiqInclude::TTest::Test(
                  x.begin(), labels.with_stride(str_labels_exp), dim_exp);
              break;
            case HetStats::MannWhitney:
              if (not_sorted) { update_sortx(); not_sorted = false; }
              pvalue_samples[z][s] = MajiqInclude::MannWhitney::Test(
                  mannwhitney, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              break;
            case HetStats::TNOM:
              if (not_sorted) { update_sortx(); not_sorted = false; }
              pvalue_samples[z][s] = MajiqInclude::TNOM::Test(
                  tnom, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              break;
            case HetStats::InfoScore:
              if (not_sorted) { update_sortx(); not_sorted = false; }
              pvalue_samples[z][s] = MajiqInclude::InfoScore::Test(
                  infoscore, x.begin(), sortx.begin(),
                  labels.with_stride(str_labels_exp), dim_exp);
              break;
            default:
              break;
          }
          // if result was nan, we know not to do it in the future
          if (std::isnan(pvalue_samples[z][s])) {
            stat_valid[z] = false;
            --ct_valid;
          }
        }  // loop to compute valid z-th statistic on psisample s
      }  // done with getting pvalues for each psisample

      // compute quantiles for stats
      for (npy_intp z = 0; z < dim_stat; ++z) {
        if (stat_valid[z]) {
          auto& in_pvalues = pvalue_samples[z];
          auto q_q = q.with_stride(str_q_q);
          for (npy_intp k = 0; k < dim_q; ++k, ++q_q) {
            using MajiqInclude::quantile;
            pvalue_quantiles[z][k] = quantile(
                in_pvalues.begin(), in_pvalues.end(), *q_q);
          }  // loop over quantiles to compute
        }  // if this stat was valid
      }  // loop over requested statistics
    }  // if any stats were valid for psisamples

    // write pvalue quantiles to output
    auto out_stat = out.with_stride(str_out_stat);
    for (npy_intp z = 0; z < dim_stat; ++z, ++out_stat) {
      if (stat_valid[z]) {
        std::copy(pvalue_quantiles[z].begin(), pvalue_quantiles[z].end(),
            out_stat.with_stride(str_out_q));
      } else {
        // invalid statistics should be NaN
        out_stat
          .with_stride(str_out_q)
          .fill(dim_q, std::numeric_limits<double>::quiet_NaN());
      }  // ifelse on whether stats were computed, set output pvalue_quantiles
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
  NPY_FLOAT, NPY_FLOAT, NPY_BOOL, NPY_FLOAT, NPY_INT64, NPY_INT64, NPY_DOUBLE, NPY_DOUBLE,
  // for use with npy_double func
  NPY_DOUBLE, NPY_DOUBLE, NPY_BOOL, NPY_DOUBLE, NPY_INT64, NPY_INT64, NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace StatsSample
}  // namespace BetaMixture
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_STATSSAMPLE_HPP
