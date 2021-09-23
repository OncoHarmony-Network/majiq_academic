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

// perform t-test on provided random-access iterators
template <typename ItX, typename ItLabels>
inline constexpr typename std::iterator_traits<ItX>::value_type
Test(ItX x, ItLabels labels, npy_intp d) {
  using RealT = typename std::iterator_traits<ItX>::value_type;
  static_assert(
      std::is_same_v<std::random_access_iterator_tag,
      typename std::iterator_traits<ItX>::iterator_category>,
      "TTest::Test requires x to be random-access iterator");
  static_assert(
      std::is_same_v<std::random_access_iterator_tag,
      typename std::iterator_traits<ItLabels>::iterator_category>,
      "TTest::Test requires labels to be random-access iterator");

  // first pass: calculate n1, n2, sum1, sum2 (ignore nan)
  int n1{0};
  int n2{0};
  RealT sum1{0};
  RealT sum2{0};
  for (npy_intp j = 0; j < d; ++j) {
    const auto& xj = x[j];
    if (npy_isnan(xj)) {
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
  for (npy_intp j = 0; j < d; ++j) {
    const auto& xj = x[j];
    if (npy_isnan(xj)) {
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
  const auto [dof, t_denom] = DOFWelchSatterthwaite(var1, var2, n1, n2);
  const RealT t = (mean1 - mean2) / t_denom;
  return TwoSidedPValue(dof, t);
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
    detail::get_value<RealT>(out_ptr) = Test(
        detail::CoreIt<RealT>::begin(x_ptr, inner_stride_x),
        detail::CoreIt<bool>::begin(labels_ptr, inner_stride_labels),
        dim_core);
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
