/*
 * Implementation of Welch's t-test (two-sample t-test assuming unequal
 * variances) for use in MAJIQ
 *
 * Copyright 2020 <University of Pennsylvania>
 *
 * Authors: Joseph K Aicher, Jorge Vaquero-Garcia
 */

#ifndef MAJIQSTATS_TTEST_HPP
#define MAJIQSTATS_TTEST_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cmath>
#include <limits>
#include <utility>

#include <boost/math/distributions/students_t.hpp>


namespace MajiqStats {

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
inline std::pair<RealT, RealT> TTestDOFWelchSatterthwaite(
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
  return std::make_pair(dof, t_denom);
}

template <typename RealT>
inline RealT TTestTwoSidedPValue(RealT dof, RealT t) {
  boost::math::students_t_distribution<RealT> t_dist(dof);
  return 2 * boost::math::cdf(t_dist, -std::fabs(t));
}

/**
 * Compute p-values for Welch's t-test on input data
 *
 * @param x: 2D array, test for each row over observations in columns
 * @param labels: 2D array, class labels for each observation
 *
 * @return p-value for test statistic for each row (1D array)
 */
template <typename RealT>
pybind11::array_t<RealT> TTest(
    pybind11::array_t<RealT> x,
    pybind11::array_t<bool> labels) {
  // check for correct number of dimensions
  if (x.ndim() != 2) {
    throw std::runtime_error("x is not 2-dimensional");
  } else if (labels.ndim() != 2) {
    throw std::runtime_error("labels is not 2-dimensional");
  }
  // check that dimensions match
  if (x.shape(0) != labels.shape(0) || x.shape(1) != labels.shape(1)) {
    throw std::runtime_error("x.shape does not match labels.shape");
  }
  //
  // create output array, 1D with length equal to rows of x
  pybind11::array_t<RealT> result(x.shape(0));
  // unchecked access to the array values
  auto _x = x.template unchecked<2>();
  auto _labels = labels.template unchecked<2>();
  auto _result = result.template mutable_unchecked<1>();

  // calculate statistics per row of x
  for (pybind11::ssize_t i = 0; i < _x.shape(0); ++i) {
    // calculate t-statistic, degrees of freedom
    int n1 = 0;
    int n2 = 0;
    RealT sum1 = RealT{0};
    RealT sum2 = RealT{0};
    // first pass: calculate n1, n2, sum1, sum2 (ignore nan)
    for (pybind11::ssize_t j = 0; j < _x.shape(1); ++j) {
      if (std::isnan(_x(i, j))) {
        // skip missing values
        continue;
      }
      if (_labels(i, j)) {
        ++n1;
        sum1 += _x(i, j);
      } else {
        ++n2;
        sum2 += _x(i, j);
      }
    }
    if (n1 < 2 || n2 < 2) {
      // not enough degrees of freedom
      _result(i) = std::numeric_limits<RealT>::quiet_NaN();
      continue;  // no need to do rest of calculation, next row
    }
    // means are then:
    const RealT mean1 = sum1 / n1;
    const RealT mean2 = sum2 / n2;
    // second pass to estimate sample variance with Bessel's correction
    RealT rss1 = RealT{0};
    RealT rss2 = RealT{0};
    for (pybind11::ssize_t j = 0; j < _x.shape(1); ++j) {
      if (std::isnan(_x(i, j))) {
        // skip missing values
        continue;
      }
      if (_labels(i, j)) {
        RealT residual = _x(i, j) - mean1;
        rss1 += residual * residual;
      } else {
        RealT residual = _x(i, j) - mean2;
        rss2 += residual * residual;
      }
    }
    const RealT var1 = rss1 / (n1 - 1);
    const RealT var2 = rss2 / (n2 - 1);
    const auto [dof, t_denom] = TTestDOFWelchSatterthwaite(var1, var2, n1, n2);
    const RealT t = (mean1 - mean2) / t_denom;
    _result(i) = TTestTwoSidedPValue(dof, t);
  }

  // return result
  return result;
}

}  // namespace MajiqStats

#endif  // MAJIQSTATS_TTEST_HPP
