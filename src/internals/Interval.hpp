/**
 * Interval.hpp
 *
 * different intervals
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_INTERVAL_HPP
#define MAJIQ_INTERVAL_HPP

#include <tuple>
#include <stdexcept>
#include <iostream>
#include <functional>

#include "MajiqTypes.hpp"


namespace majiq {
namespace detail {

struct Interval {
 public:
  position_t start;
  position_t end;

  // constructors
  Interval(position_t a, position_t b) : start{a}, end{b} { }
  Interval() : Interval{-1, -1} { }
  Interval(const Interval&) = default;
  Interval(Interval&&) = default;
  Interval& operator=(const Interval&) = default;
  Interval& operator=(Interval&&) = default;
};

// equality
inline bool operator==(const Interval& x, const Interval& y) noexcept {
  return x.start == y.start && x.end == y.end;
}
inline bool operator!=(const Interval& x, const Interval& y) noexcept {
  return !(x == y);
}

}  // namespace detail


struct OpenInterval;
struct ClosedInterval;
// closed, but treat negative coordinates as missing (i.e., for half exons)
struct ClosedOrHalfInterval;

/**
 * Open intervals
 */
struct OpenInterval : public detail::Interval {
 public:
  bool is_full_interval() const { return true; }
  // length, containment of individual coordinates
  position_t length() const { return -1 + end - start; }
  bool contains(position_t x) const { return start < x && x < end; }

  // constructors
  OpenInterval(position_t a, position_t b) : detail::Interval(a, b) {
    if (length() < 0) {
      throw std::invalid_argument("Interval must have nonnegative length");
    }
  }
  OpenInterval() : OpenInterval{-2, -1} { }
  OpenInterval(const OpenInterval&) = default;
  OpenInterval(OpenInterval&&) = default;
  OpenInterval& operator=(const OpenInterval&) = default;
  OpenInterval& operator=(OpenInterval&&) = default;
  static OpenInterval FromStartLength(position_t start, position_t length) {
    return OpenInterval{start, start + length + 1};
  }
  ClosedInterval AsClosed() const;

  // valid positions (aliases for start/end)
  const position_t& first_pos() const noexcept { return start; }
  const position_t& last_pos() const noexcept { return end; }

  // tuple representation
  std::tuple<const position_t&, const position_t&> as_tuple() const noexcept {
    return std::tie(first_pos(), last_pos());
  }
  std::tuple<const position_t&, const position_t&> rev_tuple() const noexcept {
    return std::tie(last_pos(), first_pos());
  }
};

/**
 * Closed intervals
 */
struct ClosedInterval : public detail::Interval {
 public:
  bool is_full_interval() const { return true; }
  // length, containment of individual coordinates
  position_t length() const { return 1 + end - start; }
  bool contains(position_t x) const { return start <= x && x <= end; }

  // constructors
  ClosedInterval(position_t a, position_t b) : detail::Interval(a, b) {
    if (length() < 0) {
      throw std::invalid_argument("Interval must have nonnegative length");
    }
  }
  ClosedInterval() : ClosedInterval{-1, -1} { }
  ClosedInterval(const ClosedInterval&) = default;
  ClosedInterval(ClosedInterval&&) = default;
  ClosedInterval& operator=(const ClosedInterval&) = default;
  ClosedInterval& operator=(ClosedInterval&&) = default;
  static ClosedInterval FromStartLength(position_t start, position_t length) {
    return ClosedInterval{start, start + length - 1};
  }
  OpenInterval AsOpen() const;

  // valid positions (aliases for start/end)
  const position_t& first_pos() const noexcept { return start; }
  const position_t& last_pos() const noexcept { return end; }

  // tuple representation
  std::tuple<const position_t&, const position_t&> as_tuple() const noexcept {
    return std::tie(first_pos(), last_pos());
  }
  std::tuple<const position_t&, const position_t&> rev_tuple() const noexcept {
    return std::tie(last_pos(), first_pos());
  }
};

inline ClosedInterval OpenInterval::AsClosed() const {
  return ClosedInterval{start + 1, end - 1};
}
inline OpenInterval ClosedInterval::AsOpen() const {
  return OpenInterval{start - 1, end + 1};
}

/**
 * Closed intervals, permitting missing (i.e., negative) coordinates
 */
struct ClosedOrHalfInterval : public detail::Interval {
 public:
  // full interval or otherwise?
  bool has_start() const { return start >= 0; }
  bool has_end() const { return end >= 0; }
  bool is_full_interval() const { return has_start() && has_end(); }
  bool is_half_interval() const { return has_start() != has_end(); }
  bool is_invalid() const { return !(has_start() || has_end()); }
  // for ordering of interval types when sorting [1, 1] < [1, -1] < [-1, 1]
  uint8_t interval_type_priority() const {
    return (has_start() ? 2 : 0) | (has_end() ? 1 : 0);
  }

  // length, containment of individual coordinates
  position_t length() const {
    return is_full_interval() ? 1 + end - start : 0;
  }
  bool contains(position_t x) const {
    return is_full_interval() && start <= x && x <= end;
  }

  // constructors
  ClosedOrHalfInterval(position_t a, position_t b) : detail::Interval(a, b) {
    if (is_full_interval() && length() < 0) {
      throw std::invalid_argument("Interval must have nonnegative length");
    }
  }
  ClosedOrHalfInterval() : ClosedOrHalfInterval{-1, -1} { }
  ClosedOrHalfInterval(const ClosedOrHalfInterval&) = default;
  ClosedOrHalfInterval(ClosedOrHalfInterval&&) = default;
  ClosedOrHalfInterval& operator=(const ClosedOrHalfInterval&) = default;
  ClosedOrHalfInterval& operator=(ClosedOrHalfInterval&&) = default;
  static ClosedOrHalfInterval FromStartLength(position_t start, position_t length) {
    return ClosedOrHalfInterval{start, start + length - 1};
  }

  // valid positions (if any)
  const position_t& first_pos() const noexcept {
    return has_start() ? start : end;
  }
  const position_t& last_pos() const noexcept {
    return has_end() ? end : start;
  }

  // tuple representation
  std::tuple<const position_t&, const position_t&> as_tuple() const noexcept {
    return std::tie(first_pos(), last_pos());
  }
  std::tuple<const position_t&, const position_t&> rev_tuple() const noexcept {
    return std::tie(last_pos(), first_pos());
  }
};

// ordering
template <
  typename I1,
  typename I2,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I2>, bool> = true
>
inline bool operator<(const I1& x, const I2& y) noexcept {
  return x.as_tuple() < y.as_tuple();
}
/**
 * explicit specialization of operator< for ClosedOrHalfInterval vs self
 *
 * In order to have a full ordering, need ordering on missing status
 */
template <>
inline bool operator<(
    const ClosedOrHalfInterval& x, const ClosedOrHalfInterval& y) noexcept {
  const auto x_priority = x.interval_type_priority();
  const auto y_priority = y.interval_type_priority();
  // x and y priority are switched since we want higher priority first
  return std::tie(x.first_pos(), x.last_pos(), y_priority)
    < std::tie(y.first_pos(), y.last_pos(), x_priority);
}
// ordering vs coordinate is against first position
template <
  typename I1,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true
>
inline bool operator<(const I1& x, const position_t& y) noexcept {
  return x.first_pos() < y;
}
template <
  typename I1,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true
>
inline bool operator<(const position_t& x, const I1& y) noexcept {
  return x < y.first_pos();
}
// derived orderings
template <
  typename I1,
  typename I2,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I2>, bool> = true
>
inline bool operator>(const I1& x, const I2& y) noexcept {
  return y < x;
}
template <
  typename I1,
  typename I2,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I2>, bool> = true
>
inline bool operator<=(const I1& x, const I2& y) noexcept {
  return !(x > y);
}
template <
  typename I1,
  typename I2,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I2>, bool> = true
>
inline bool operator>=(const I1& x, const I2& y) noexcept {
  return !(x < y);
}

// how to print intervals
template <
  typename I1,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true
>
inline std::ostream& operator<<(std::ostream& os, const I1& x) {
  os << x.start << "-" << x.end;
  return os;
}
template <>
inline std::ostream& operator<<(
    std::ostream& os, const ClosedOrHalfInterval& x) {
  if (x.has_start()) { os << x.start; } else { os << "na"; }
  os << "-";
  if (x.has_end()) { os << x.end; } else { os << "na"; }
  return os;
}

// subset/superset of intervals
template <
  typename I1,
  typename I2,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I2>, bool> = true
>
inline bool IntervalSubsets(const I1& sub, const I2& sup) noexcept {
  if constexpr(std::is_base_of_v<ClosedOrHalfInterval, I1>) {
    if (!sub.is_full_interval()) { return false; }
  }
  return sup.contains(sub.start) && sup.contains(sub.end);
}

// interval preceding
// NOTE: we consider closed-open to be intersecting if they share a coordinate.
// That is, (2, 4) intersects [4, 6]. This reflects our use of open intervals
// for junctions and closed intervals for exons.
template <
  typename I1,
  typename I2,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I2>, bool> = true
>
inline bool IntervalPrecedes(const I1& before, const I2& after) noexcept {
  const auto& before_last = before.last_pos();
  return before_last < after.first_pos() && before_last < after.last_pos();
}
// if both are open, then values can be equal and still be 'preceding'
template <>
inline bool IntervalPrecedes(
    const OpenInterval& before, const OpenInterval& after) noexcept {
  return before.end <= after.start;
}

struct IntervalPrecedesT {
  template <typename T, typename U>
  inline bool operator()(const T& before, const U& after) const noexcept {
    return IntervalPrecedes<T, U>(before, after);
  }
};

// intersection of intervals
template <typename I1, typename I2>
inline bool IntervalIntersects(const I1& x, const I2& y) noexcept {
  return !(IntervalPrecedes(x, y) || IntervalPrecedes(y, x));
}

}  // namespace majiq


#endif  // MAJIQ_INTERVAL_HPP
