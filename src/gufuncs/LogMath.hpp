/**
 * LogMath.hpp
 *
 * Helper functions for doing combinatorics/addition in logspace
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQGUFUNCS_LOGMATH_HPP
#define MAJIQGUFUNCS_LOGMATH_HPP

#include <cmath>
#include <limits>


namespace MajiqGufuncs {
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

}  // namespace detail
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_LOGMATH_HPP
