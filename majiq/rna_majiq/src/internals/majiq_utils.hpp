/**
 * majiq_utils.hpp
 *
 * Helper functions for MAJIQ
 *
 * Author: Joseph K Aicher
 */
#ifndef MAJIQ_UTILS_HPP
#define MAJIQ_UTILS_HPP

#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>
#include <vector>

namespace majiq {

/**
 * Obtains median of values between first and last
 *
 * @param first, last random-access iterators. Note: values will be reordered
 */
template<bool sorted = false, class It>
inline
#if __cplusplus >= 201703L
constexpr  // required for c++17 and later
#endif
typename std::iterator_traits<It>::value_type median(It first, It last) {
  using T = typename std::iterator_traits<It>::value_type;
  // should be random iterator
  static_assert(std::is_same<
      std::random_access_iterator_tag,
      typename std::iterator_traits<It>::iterator_category>::value,
      "median() only accepts random-access iterators as input\n");
  // how many elements?
  const auto n = std::distance(first, last);
  // 3 cases
  if (n < 1) {
    return std::numeric_limits<T>::quiet_NaN();
  } else if (n == 1) {
    return *first;
  } else if (n == 2) {
      // no need to sort anything
    return (first[0] + first[1]) / 2;
  } else {
    // sort to get the value that would be in the middle if sorted
    It it_below = first + ((n - 1) / 2);
    // make sure that value for it_below avaiable as if sorted
#if __cplusplus >= 201703L
    if constexpr(!sorted) {
#else
    if (!sorted) {
#endif
      std::nth_element(first, it_below, last);
    }
    if (n % 2 == 1) {
      // odd means that this is the median
      return *it_below;
    } else {
      // we need the value that follows in a sorted array and take the average
#if __cplusplus >= 201703L
      if constexpr(sorted) {
#else
      if (sorted) {
#endif
        // sorted --> it_below and it_below + 1 point to what to use
        return (it_below[0] + it_below[1]) / 2;
      } else {
        // next element is the minimum value after it_below
        return (*it_below + *std::min_element(it_below + 1, last)) / 2;
      }
    }
  }
}

/**
 * Obtains the specified quantile of values between first and last
 *
 * @param first, last random-access iterators. Note: values will be reordered
 * @param q the quantile to evaluate
 *
 * Obtains quantile between observations using linear interpolation (i.e. NumPy
 * default)
 */
template<bool sorted = false, class It>
inline
#if __cplusplus >= 201703L
constexpr  // required for c++17 and later
#endif
typename std::iterator_traits<It>::value_type quantile(
    It first, It last, float q) {
  using T = typename std::iterator_traits<It>::value_type;
  // should be random iterator
  static_assert(std::is_same<
      std::random_access_iterator_tag,
      typename std::iterator_traits<It>::iterator_category>::value,
      "quantile() only accepts random-access iterators as input\n");
  // how many elements?
  const auto n = std::distance(first, last);
  if (n < 1) {
    return std::numeric_limits<T>::quiet_NaN();
  } else if (n == 1) {
    return *first;
  } else {
    // otherwise, we need to determine where we belong between 0 and n-1
    if (q <= 0) {
#if __cplusplus >= 201703L
      if constexpr(sorted) {
#else
      if (sorted) {
#endif
        return *first;
      } else {
        return *std::min_element(first, last);
      }
    } else if (q >= 1) {
#if __cplusplus >= 201703L
      if constexpr(sorted) {
#else
      if (sorted) {
#endif
        return *(last - 1);
      } else {
        return *std::max_element(first, last);
      }
    } else {
      // between which indexes in sorted array is our quantile?
      const float idx_float = (n - 1) * q;
      const float idx_floor_float = std::floor(idx_float);
      using dIt = typename std::iterator_traits<It>::difference_type;
      It it_floor = first + static_cast<dIt>(idx_floor_float);
      // declare values for interpolation used in quantile
      T value_below, value_above;
#if __cplusplus >= 201703L
      if constexpr(sorted) {
#else
      if (sorted) {
#endif
        value_below = it_floor[0];
        value_above = it_floor[1];
      } else {
        // partial sort of values to put correct value at it_floor
        std::nth_element(first, it_floor, last);
        value_below = *it_floor;
        // value above is the smallest subsequent element
        value_above = *std::min_element(it_floor + 1, last);
      }
      // linear interpolation between the two:
      return value_below
        + (idx_float - idx_floor_float) * (value_above - value_below);
    }
  }
}

/**
 * Obtains mean of values between first and last
 *
 * @param first, last iterators
 */
template<class It>
inline
#if __cplusplus >= 201703L
constexpr  // required for c++17 and later
#endif
typename std::iterator_traits<It>::value_type mean(It first, It last) {
  using T = typename std::iterator_traits<It>::value_type;
  const auto n = std::distance(first, last);
  if (n < 1) {
    return std::numeric_limits<T>::quiet_NaN();
  } else {
    return std::accumulate(first, last, T{}) / n;
  }
}

/**
 * Obtains variance of values between first and last
 *
 * @param first, last iterators
 * @param known_mean mean of data between first and last
 * @param ddof template parameter for Bessel's correction
 */
template<int ddof = 1, class It>
inline
#if __cplusplus >= 201703L
constexpr  // required for c++17 and later
#endif
typename std::iterator_traits<It>::value_type variance(
    It first, It last,
    typename std::iterator_traits<It>::value_type known_mean) {
  using T = typename std::iterator_traits<It>::value_type;
  const auto n = std::distance(first, last);
  if (n <= ddof) {
    return std::numeric_limits<T>::quiet_NaN();
  } else {
    // get residual sum of squares, divide by n adjusted by ddof
    return std::accumulate(first, last, T{}, [known_mean](T rss, T x) {
          T residual = x - known_mean;
          return rss + residual * residual;
        }) / (n - ddof);
  }
}

/**
 * variance of input values without knowledge of mean
 */
template<int ddof = 1, class It>
inline
#if __cplusplus >= 201703L
constexpr  // required for c++17 and later
#endif
typename std::iterator_traits<It>::value_type variance(It first, It last) {
  using T = typename std::iterator_traits<It>::value_type;
  const T known_mean = mean(first, last);
  return variance<ddof, It>(first, last, known_mean);
}

/**
 * logsumexp of input values
 */
template<class It>
inline
#if __cplusplus >= 201703L
constexpr  // required for c++17 and later
#endif
typename std::iterator_traits<It>::value_type logsumexp(It first, It last) {
  using T = typename std::iterator_traits<It>::value_type;
  const auto n = std::distance(first, last);
  if (n < 1) {
    return std::numeric_limits<T>::quiet_NaN();
  } else if (n == 1) {
    return *first;
  } else {
    // get maximum value
    const T max_elem = *std::max_element(first, last);
    // get sum of exp(x - max_elem)
    const T sum_exp = std::accumulate(first, last, T{},
        [max_elem](T sum, T x) {
          return sum + std::exp(x - max_elem);
        });
    // return log of sum_exp, but add back max_elem to result
    return std::log(sum_exp) + max_elem;
  }
}

/**
 * logsumexp over vector of vectors (2D)
 *
 * @note a 2D matrix really should be stored contiguously (vector of vectors
 * is *not* contiguous), in which case we could just use logsumexp. But
 * implemented for historical reasons
 */
template<typename T>
inline T logsumexp_2D(std::vector<std::vector<T>> values) {
  if (values.empty()) {
    return std::numeric_limits<T>::quiet_NaN();
  } else if (values.size() == 1) {
    // only one vector<T> to operate over
    return logsumexp(values[0].begin(), values[0].end());
  } else {
    // get maximum value (initialize with max from first vector)
    T max_elem = *std::max_element(values[0].begin(), values[0].end());
    max_elem = std::accumulate(values.begin() + 1, values.end(), max_elem,
        [](T cur_max, std::vector<T> v) -> T {
          const T v_max = *std::max_element(v.begin(), v.end());
          return std::max(cur_max, v_max);
        });
    // get sum of exp(x - max_elem)
    const T sum_exp = std::accumulate(values.begin(), values.end(), T{},
        [max_elem](T outer_sum, std::vector<T> v) -> T {
          return outer_sum + std::accumulate(v.begin(), v.end(), T{},
              [max_elem](T inner_sum, T x) -> T {
                return inner_sum + std::exp(x - max_elem);
              });
        });
    // return log of sum_exp but add back max_elem
    return std::log(sum_exp) + max_elem;
  }
}

}  // namespace majiq

#endif  // MAJIQ_UTILS_HPP
