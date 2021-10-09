/**
 * TTest.hpp
 *
 * Implementation of Welch's t-test (two-sample t-test assuming unequal
 * variances) for use in MAJIQ
 *
 * Copyright 2020 <University of Pennsylvania>
 *
 * Authors: Joseph K Aicher, Jorge Vaquero-Garcia
 */

#ifndef MAJIQINCLUDE_STATS_TTEST_HPP
#define MAJIQINCLUDE_STATS_TTEST_HPP

#include <cmath>
#include <limits>
#include <utility>

#include <boost/math/distributions/students_t.hpp>


namespace MajiqInclude {
namespace TTest {

template <typename RealT>
struct DOF_pair {
  RealT dof;
  RealT t_denom;
};

/** Estimate pooled degrees of freedom for Welch t-test
 *
 * Estimate pooled degrees of freedom for Welch t-test using
 * Welch-Satterthwaite equation
 *
 * @param var1 sample variance of first sample
 * @param var2 sample variance of second sample
 * @param n1 sample size of first sample
 * @param n2 sample size of second sample
 *
 * @note Modeled after _unequal_var_ttest_denom in SciPy
 *
 * @returns degrees of freedom, denominator for t-statistic
 */
template <typename RealT>
inline DOF_pair<RealT> DOFWelchSatterthwaite(
    RealT var1, RealT var2, int n1, int n2) {
  // compute ratio of variance to sample size
  const RealT ratio1 = var1 / n1;
  const RealT ratio2 = var2 / n2;
  // compute components of dof
  const RealT numerator_term = ratio1 + ratio2;
  const RealT denom1 = ratio1 * ratio1 / (n1 - 1);
  const RealT denom2 = ratio2 * ratio2 / (n2 - 1);
  // estimated degrees of freedom are
  const RealT dof = numerator_term * numerator_term / (denom1 + denom2);
  // denominator for t-test is
  const RealT t_denom = std::sqrt(ratio1 + ratio2);
  // return pair of dof, t_denom
  return DOF_pair<RealT>{dof, t_denom};
}

template <typename RealT>
inline RealT TwoSidedPValue(RealT dof, RealT t) {
  boost::math::students_t_distribution<RealT> t_dist(dof);
  return 2 * boost::math::cdf(t_dist, -std::fabs(t));
}

// perform t-test on provided random-access iterators
template <typename ItX, typename ItLabels,
         typename RealT = typename std::iterator_traits<ItX>::value_type>
inline RealT Test(ItX x, ItLabels labels, int64_t d) {
  // first pass: calculate n1, n2, sum1, sum2 (ignore nan)
  int n1{0};
  int n2{0};
  RealT sum1{0};
  RealT sum2{0};
  for (int64_t j = 0; j < d; ++j) {
    const auto& xj = x[j];
    if (std::isnan(xj)) {
      // skip missing values
      continue;
    }
    if (labels[j]) {
      ++n1;
      sum1 += xj;
    } else {
      ++n2;
      sum2 += xj;
    }
  }  // first pass over core dimension for sum1/2
  if (n1 < 2 || n2 < 2) {
    // not enough degrees of freedom, return NaN
    return std::numeric_limits<RealT>::quiet_NaN();
  }
  const RealT mean1 = sum1 / n1;
  const RealT mean2 = sum2 / n2;

  // second pass to estimate sample variance with Bessel's correction
  RealT rss1{0};
  RealT rss2{0};
  for (int64_t j = 0; j < d; ++j) {
    const auto& xj = x[j];
    if (std::isnan(xj)) {
      // skip missing values
      continue;
    }
    if (labels[j]) {
      const RealT residual = xj - mean1;
      rss1 += residual * residual;
    } else {
      const RealT residual = xj - mean2;
      rss2 += residual * residual;
    }
  }  // second pass over core dimension for rss1/2
  const RealT var1 = rss1 / (n1 - 1);
  const RealT var2 = rss2 / (n2 - 1);
  const auto dof_pair = DOFWelchSatterthwaite(var1, var2, n1, n2);
  const RealT t = (mean1 - mean2) / dof_pair.t_denom;
  return TwoSidedPValue(dof_pair.dof, t);
}

}  // namespace TTest
}  // namespace MajiqInclude

#endif  // MAJIQINCLUDE_STATS_TTEST_HPP
