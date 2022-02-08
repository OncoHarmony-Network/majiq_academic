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

/*
 * Determine if interval `before` precedes `after`.
 *
 * If neither precedes the other, the intervals overlap/intersect.
 *
 * Special care needs to be taken for the special case of zero-length
 * intervals.
 * When defining introns, we would like to be able to have directly adjacent
 * zero length introns: for example, [100, 99], [100, 200], [201, 200].
 * The zero length introns in this case would correspond to unspliced reads
 * over the adjoining exon boundaries.
 * For this case, we consider zero length regions as non-overlapping.
 * However, when assigning coverage to introns, or matching introns, we would
 * like these introns to be considered as overlapping.
 * For example, if the intron [100, 200] is split by a denovo exon [100, 200],
 * we would get two zero-length introns [100, 99] and [201, 200] that need to
 * be matched back to the original intron [100, 200].
 * For this case, the zero length regions are overlapping.
 *
 * We use the boolean template parameter bool ZERO_REGION_OVERLAPS to handle
 * this case.
 *
 * We consider closed-open to be intersecting if they share a coordinate. That
 * is, (2, 4) intersects [4, 6], reflecting our use of open intervals for
 * junctions and closed intervals for exons.
 */
template <
  bool ZERO_REGION_OVERLAPS = true,
  typename I1,
  typename I2,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I1>, bool> = true,
  std::enable_if_t<std::is_base_of_v<detail::Interval, I2>, bool> = true
>
inline bool IntervalPrecedes(const I1& before, const I2& after) noexcept {
  if constexpr(std::is_same_v<I1, OpenInterval>
      && std::is_same_v<I2, OpenInterval>) {
    // treat open intervals as special case to allow overlapping end/start
    return before.end <= after.start;
  } else {
    // length() treats half intervals as zero length, so we use this condition
    // this condition isn't semantically correct for open intervals, but it is
    // correct for this function, especially since we specialize for
    // OpenInterval vs OpenInterval above
    const auto& before_first = before.first_pos();
    const auto& before_last = before.last_pos();
    const auto& after_first = after.first_pos();
    const auto& after_last = after.last_pos();
    const bool before_zero = before_first > before_last;
    const bool after_zero = after_first > after_last;
    if (before_zero && after_zero) {
      // if both are zero length, then should check start vs start or end vs end
      return before_last < after_last;
    } else if (before_zero) {
      // - after is not zero length (would have been first case)
      // - before_last < before_first
      if constexpr(ZERO_REGION_OVERLAPS) {
        // [a, a - 1] DOES NOT precede [a, b] (b >= a)
        return before_first < after_first;
      } else {
        // [a, a - 1] DOES precede [a, b] (b >= a)
        return before_first <= after_first;
      }
    } else if (after_zero) {
      // - before is not zero length (would have been first case)
      // - after_last < after_first
      if constexpr(ZERO_REGION_OVERLAPS) {
        // [a, b] (a <= b) DOES NOT precede [b + 1, b]
        return before_last < after_last;
      } else {
        // [a, b] (a <= b) DOES precede [b + 1, b]
        return before_last <= after_last;
      }
    } else {
      // both have length and we can compare last vs first as expected
      return before_last < after_first;
    }
  }
}

template <bool ZERO_REGION_OVERLAPS = true>
struct IntervalPrecedesT {
  template <typename T, typename U>
  inline bool operator()(const T& before, const U& after) const noexcept {
    return IntervalPrecedes<ZERO_REGION_OVERLAPS>(before, after);
  }
};

// intersection of intervals
template <bool ZERO_REGION_OVERLAPS = true, typename I1, typename I2>
inline bool IntervalIntersects(const I1& x, const I2& y) noexcept {
  return !(IntervalPrecedes<ZERO_REGION_OVERLAPS>(x, y)
      || IntervalPrecedes<ZERO_REGION_OVERLAPS>(y, x));
}

}  // namespace majiq


#endif  // MAJIQ_INTERVAL_HPP
