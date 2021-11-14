/**
 * histogram.hpp
 *
 * Implementation of histogram with uniform bins
 *
 * Copyright 2021 <University of Pennsylvania>
 */

#ifndef MAJIQINCLUDE_HISTOGRAM_HPP
#define MAJIQINCLUDE_HISTOGRAM_HPP

#include <algorithm>
#include <cmath>

namespace MajiqInclude {
namespace Histogram {

template <typename ItX, typename ItOut, typename RealT,
         typename CountT = typename std::iterator_traits<ItOut>::value_type,
         typename std::enable_if<
          std::is_integral<CountT>::value, bool>::type = true,
         typename std::enable_if<
          std::is_floating_point<RealT>::value, bool>::type = true,
         typename std::enable_if<
          std::is_same<
            CountT, typename std::iterator_traits<ItOut>::value_type
          >::value
         , bool>::type = true,
         typename std::enable_if<
          std::is_same<
            RealT, typename std::iterator_traits<ItX>::value_type
          >::value
         , bool>::type = true>
inline void histogram(
    ItX x, const int64_t n_x, RealT min_x, RealT max_x,
    ItOut out, const int64_t n_out) {
  // fill out with zeros
  if (n_out < 1) {
    return;
  } else if (n_out > 1) {
    std::fill(out, out + n_out, CountT{0});
  } else {
    *out = 0;
  }
  // if invalid range, no need to count, each bin will be zero
  if (std::isnan(min_x) || std::isnan(max_x) || min_x >= max_x) {
    return;
  }
  // count values of x that fall within range [min_x, max_x)
  for (int64_t i = 0; i < n_x; ++i, ++x) {
    if (std::isnan(*x)) { continue; }
    auto index = static_cast<int64_t>(std::floor(
                n_out * (*x - min_x) / (max_x - min_x)));
    if (index < 0 || index >= n_out) { continue; }
    ++out[index];
  }
  return;
}

}  // namespace Histogram
}  // namespace MajiqInclude


#endif  // MAJIQINCLUDE_HISTOGRAM_HPP
