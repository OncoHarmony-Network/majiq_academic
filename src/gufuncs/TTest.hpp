/*
 * Implementation of Welch's t-test (two-sample t-test assuming unequal
 * variances) for use in MAJIQ
 *
 * Copyright 2020 <University of Pennsylvania>
 *
 * Authors: Joseph K Aicher, Jorge Vaquero-Garcia
 */

#ifndef MAJIQGUFUNCS_TTEST_HPP
#define MAJIQGUFUNCS_TTEST_HPP

#include <numpy/ndarraytypes.h>

#include <cmath>
#include <limits>
#include <utility>

#include "helpers.hpp"

#include <boost/math/distributions/students_t.hpp>


namespace MajiqGufuncs {
namespace TTest {

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
inline std::pair<RealT, RealT> DOFWelchSatterthwaite(
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
inline RealT TwoSidedPValue(RealT dof, RealT t) {
  boost::math::students_t_distribution<RealT> t_dist(dof);
  return 2 * boost::math::cdf(t_dist, -std::fabs(t));
}

template <typename RealT>
inline void Inner(
    char* x, char* labels, char* out, const npy_intp d,
    const npy_intp s_x, const npy_intp s_labels) {
  using detail::get_value;
  // calculate t-statistic, degrees of freedom
  int n1{0};
  int n2{0};
  RealT sum1{0};
  RealT sum2{0};
  // first pass: calculate n1, n2, sum1, sum2 (ignore nan)
  for (npy_intp j = 0; j < d; ++j) {
    const RealT& xi = get_value<RealT>(x, j, s_x);
    if (npy_isnan(xi)) {
      // skip missing values
      continue;
    }
    const bool& labelsi = get_value<bool>(labels, j, s_labels);
    if (labelsi) {
      ++n1;
      sum1 += xi;
    } else {
      ++n2;
      sum2 += xi;
    }
  }  // first pass over core dimension for sum1/2
  if (n1 < 2 || n2 < 2) {
    // not enough degrees of freedom, return NaN
    get_value<RealT>(out, 0, 0) = std::numeric_limits<RealT>::quiet_NaN();
    return;
  }
  const RealT mean1 = sum1 / n1;
  const RealT mean2 = sum2 / n2;
  // second pass to estimate sample variance with Bessel's correction
  RealT rss1{0};
  RealT rss2{0};
  for (npy_intp j = 0; j < d; ++j) {
    const RealT& xi = get_value<RealT>(x, j, s_x);
    if (npy_isnan(xi)) {
      // skip missing values
      continue;
    }
    const bool& labelsi = get_value<bool>(labels, j, s_labels);
    if (labelsi) {
      const RealT residual = xi - mean1;
      rss1 += residual * residual;
    } else {
      const RealT residual = xi - mean2;
      rss2 += residual * residual;
    }
  }  // second pass over core dimension for rss1/2
  const RealT var1 = rss1 / (n1 - 1);
  const RealT var2 = rss2 / (n2 - 1);
  const auto [dof, t_denom] = DOFWelchSatterthwaite(var1, var2, n1, n2);
  const RealT t = (mean1 - mean2) / t_denom;
  get_value<RealT>(out, 0, 0) = TwoSidedPValue(dof, t);
  return;
}

template <typename RealT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides on each variable for outer loop
  const npy_intp stride_x = *steps++;
  const npy_intp stride_labels = *steps++;
  const npy_intp stride_out = *steps++;
  // core dimensions
  const npy_intp dim_core = dimensions[0];
  // inner strides
  const npy_intp inner_stride_x = steps[0];
  const npy_intp inner_stride_labels = steps[1];
  // pointers to data
  char* x_ptr = args[0];
  char* labels_ptr = args[1];
  char* out_ptr = args[2];

  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i,
      x_ptr += stride_x, labels_ptr += stride_labels, out_ptr += stride_out) {
    Inner<RealT>(x_ptr, labels_ptr, out_ptr,
        dim_core, inner_stride_x, inner_stride_labels);
  }
  return;
}

static char name[] = "ttest";
static char signature[] = "(n),(n)->()";
static char doc[] = R"pbdoc(
Compute p-values for Welch's t-test on input data

Compute p-values for Welch's t-test on input data, using a two-sided
alternative hypothesis and omitting nan values.

Parameters
----------
x: array[float]
    test over observations in last axis
labels: array[bool]
    test over labels in last axis

Returns
-------
array[float]
    broadcast p-values for observations/labels. Invalid tests are nan.
)pbdoc";
constexpr int ntypes = 2;
constexpr int nin = 2;
constexpr int nout = 1;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>)
};
static char types[
  ntypes * (nin + nout)
] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_BOOL, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_BOOL, NPY_DOUBLE,
};


}  // namespace TTest
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_TTEST_HPP
