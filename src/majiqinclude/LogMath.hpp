/**
 * LogMath.hpp
 *
 * Helper functions for doing combinatorics/addition in logspace
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQINCLUDE_LOGMATH_HPP
#define MAJIQINCLUDE_LOGMATH_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>


namespace MajiqInclude {
namespace detail {

/**
 * calculate log number of combinations: log(N choose k)
 *
 * @param k: number chosen
 * @param N: population chosen from
 */
inline double lchoose(int k, int N) {
  return std::lgamma(N + 1) - (std::lgamma(k + 1) + std::lgamma(N - k + 1));
}

/**
 * add two numbers in logspace together
 *
 * @param logx, logy
 *
 * @return log(x + y)
 */
template <typename RealT>
inline RealT logadd(RealT logx, RealT logy) {
  // if it is as small as it can get, return the other number
  if (logx <= std::numeric_limits<RealT>::lowest()) {
    return logy;
  } else if (logy <= std::numeric_limits<RealT>::lowest()) {
    return logx;
  }
  // otherwise
  if (logx >= logy) {
    logy -= logx;  // smaller number relative to larger number
  } else {
    RealT t = logx;  // the smaller number
    logx = logy;  // make logx the larger number
    logy = t - logx;  // smaller number relative to larger number
  }
  return logx + std::log1p(std::exp(logy));
}

/**
 * logsumexp of input values
 */
template <typename ItX,
         typename RealT = typename std::iterator_traits<ItX>::value_type>
inline RealT logsumexp(ItX values, const int64_t n) {
  if (n < 1) {
    return std::numeric_limits<RealT>::quiet_NaN();
  } else if (n == 1) {
    return *values;
  } else {
    auto last = values + n;
    const RealT max_val = *std::max_element(values, last);
    const RealT sum_exp = std::accumulate(values, last, RealT{},
        [max_val](RealT sum, RealT x) { return sum + std::exp(x - max_val); });
    return std::log(sum_exp) + max_val;
  }
}

/**
 * normalize log probabilities
 */
template <typename ItP>
inline void lognormalize(ItP logp, const int64_t n) {
  const auto logZ = logsumexp(logp, n);
  for (int64_t i = 0; i < n; ++i, ++logp) {
    *logp -= logZ;
  }
  return;
}

/**
 * Compute logsumexp on diagonal of A[i, j] = values_x[i] + values_y[j] without
 * computing A[i, j]
 */
template <typename ItX, typename ItY,
         typename RealT = typename std::iterator_traits<ItX>::value_type>
inline RealT logsumexp_diag(ItX values_x, ItY values_y,
    const int64_t n, const int64_t offset) {
  if (offset < 0) {
    values_x -= offset;
  } else {
    values_y += offset;
  }
  const int64_t diagonal_n = n - (offset > 0 ? offset : -offset);
  if (diagonal_n < 1) {
    return std::numeric_limits<RealT>::quiet_NaN();
  } else if (diagonal_n == 1) {
    return *values_x + *values_y;
  } else {
    RealT max_val = std::numeric_limits<RealT>::lowest();
    {  // get maximum value on diagonal
      for (int64_t i = 0; i < diagonal_n; ++i) {
        max_val = std::max(max_val, values_x[i] + values_y[i]);
      }
    }
    RealT sum_val{0.};
    for (int64_t i = 0; i < diagonal_n; ++i) {
      sum_val += std::exp(values_x[i] + values_y[i] - max_val);
    }
    return std::log(sum_val) + max_val;
  }
}

}  // namespace detail
}  // namespace MajiqInclude

#endif  // MAJIQINCLUDE_LOGMATH_HPP
