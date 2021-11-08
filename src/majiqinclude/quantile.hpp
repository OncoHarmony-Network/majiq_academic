/**
 * quantile.hpp
 *
 * Implementation of quantile on random-access iterators
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQINCLUDE_QUANTILE_HPP
#define MAJIQINCLUDE_QUANTILE_HPP

#include <algorithm>
#include <cmath>
#include <limits>

namespace MajiqInclude {

/**
 * Compute quantiles from unsorted random-access iterator. Reorders input.
 * Linear interpolation between values.
 */
template <typename It,
         typename T = typename std::iterator_traits<It>::value_type>
inline T quantile(It first, It last, float q) {
  const auto n = std::distance(first, last);
  if (n < 1) {
    return std::numeric_limits<T>::quiet_NaN();
  } else if (n == 1) {
    return *first;
  } else {
    if (q <= 0) {
      return *std::min_element(first, last);
    } else if (q >= 1) {
      return *std::max_element(first, last);
    } else {
      // between which indexes in sorted array is our quantile?
      const float idx_float = (n - 1) * q;
      const float idx_floor_float = std::floor(idx_float);
      using dIt = typename std::iterator_traits<It>::difference_type;
      It it_floor = first + static_cast<dIt>(idx_floor_float);
      // partial sort of values to put correct value at it_floor
      std::nth_element(first, it_floor, last);
      T value_below = *it_floor;
      T value_above = *std::min_element(it_floor + 1, last);
      return value_below
        + (idx_float - idx_floor_float) * (value_above - value_below);
    }
  }
}

/**
 * Median as specific case of quantile
 */
template <typename It,
         typename T = typename std::iterator_traits<It>::value_type>
inline T median(It first, It last) {
  const auto n = std::distance(first, last);
  if (n < 1) {
    return std::numeric_limits<T>::quiet_NaN();
  } else if (n == 1) {
    return *first;
  } else if (n == 2) {
    return (*first + *(++first)) / T{2};
  } else {
    // midpoint (n odd), or one more than midpoint (n even)
    const auto idx_midpoint = n / 2;
    It midpoint = first + idx_midpoint;
    std::nth_element(first, midpoint, last);
    if (n % 2 == 1) {  // odd, this is our median
      return *midpoint;
    } else {
      // actually half index past where midpoint would be
      return (*midpoint + *std::max_element(first, midpoint)) / T{2};
    }
  }
}

}  // namespace MajiqInclude

#endif  // MAJIQINCLUDE_QUANTILE_HPP
