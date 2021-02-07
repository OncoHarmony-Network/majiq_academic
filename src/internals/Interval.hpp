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
#include <boost/functional/hash.hpp>

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

  // full interval or otherwise?
  bool has_start() const { return start >= 0; }
  bool has_end() const { return end >= 0; }
  bool is_full_interval() const { return has_start() && has_end(); }
  bool is_half_interval() const { return has_start() != has_end(); }
  bool is_invalid() const { return !(has_start() || has_end()); }

  // valid positions (if any)
  const position_t& first_pos() const noexcept {
    return has_start() ? start : end;
  }
  const position_t& last_pos() const noexcept {
    return has_end() ? end : start;
  }
};
// ordering
inline bool operator<(const Interval& x, const Interval& y) noexcept {
  return std::tie(x.first_pos(), x.last_pos())
    < std::tie(y.first_pos(), y.last_pos());
}
// ordering vs coordinate is against first position
inline bool operator<(const Interval& x, const position_t& y) noexcept {
  return x.first_pos() < y;
}
inline bool operator<(const position_t& x, const Interval& y) noexcept {
  return x < y.first_pos();
}
// derived orderings
inline bool operator>(const Interval& x, const Interval& y) noexcept {
  return y < x;
}
inline bool operator<=(const Interval& x, const Interval& y) noexcept {
  return !(x > y);
}
inline bool operator>=(const Interval& x, const Interval& y) noexcept {
  return !(x < y);
}
// equality
inline bool operator==(const Interval& x, const Interval& y) noexcept {
  return x.start == y.start && x.end == y.end;
}
inline bool operator!=(const Interval& x, const Interval& y) noexcept {
  return !(x == y);
}
}  // namespace detail


/**
 * Closed intervals
 */
struct ClosedInterval : public detail::Interval {
 public:
  // length, containment of individual coordinates
  position_t length() const {
    return is_full_interval() ? 1 + end - start : 0;
  }
  bool contains(position_t x) {
    return is_full_interval() && start <= x && x <= end;
  }

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
};

/**
 * Open intervals
 */
struct OpenInterval : public detail::Interval {
 public:
  // length, containment of individual coordinates
  position_t length() const {
    return is_full_interval() ? -1 + end - start : 0;
  }
  bool contains(position_t x) {
    return is_full_interval() && start < x && x < end;
  }

  // constructors
  OpenInterval(position_t a, position_t b) : detail::Interval(a, b) {
    if (length() < 0) {
      throw std::invalid_argument("Interval must have nonnegative length");
    }
  }
  OpenInterval() : OpenInterval{-1, -1} { }
  OpenInterval(const OpenInterval&) = default;
  OpenInterval(OpenInterval&&) = default;
  OpenInterval& operator=(const OpenInterval&) = default;
  OpenInterval& operator=(OpenInterval&&) = default;
};

// intersection of intervals
template <class I1, class I2>
inline bool IntervalIntersects(const I1& x, const I2& y) {
  return x.is_full_interval() && y.is_full_interval()
    && x.start < y.end && y.start < x.end;
}
template <>  // explicit override for closed vs closed
inline bool IntervalIntersects(
    const ClosedInterval& x, const ClosedInterval& y) {
  return !(x.is_invalid() || y.is_invalid())
    && x.first_pos() <= y.last_pos() && y.first_pos() <= x.last_pos();
}

// subset/superset of intervals when they are the same type
template <class IntervalT>
inline bool IntervalSubsets(const IntervalT& sub, const IntervalT& sup) {
  return sub.is_full_interval() && sup.is_full_interval()
    && sup.start <= sub.start && sub.end <= sup.end;
}
inline bool IntervalSubsets(
    const OpenInterval& sub, const ClosedInterval& sup) {
  return sub.is_full_interval() && sup.is_full_interval()
    // adjust for open interval
    && sup.start <= sub.start + 1 && sub.end - 1 <= sup.end;
}
inline bool IntervalSubsets(
    const ClosedInterval& sub, const OpenInterval& sup) {
  return sub.is_full_interval() && sup.is_full_interval()
    // adjust for open interval
    && sup.start + 1 <= sub.start && sub.end <= sup.end - 1;
}

// how to print intervals
inline std::ostream& operator<<(std::ostream& os, const ClosedInterval& x) {
  os << "[" << x.start << ", " << x.end << "]";
  return os;
}
inline std::ostream& operator<<(std::ostream& os, const OpenInterval& x) {
  os << "(" << x.start << ", " << x.end << ")";
  return os;
}

// how to hash intervals
inline std::size_t hash_value(const OpenInterval& x) noexcept {
  std::size_t result = boost::hash_value(x.start);
  boost::hash_combine(result, x.end);
  return result;
}
inline std::size_t hash_value(const ClosedInterval& x) noexcept {
  std::size_t result = boost::hash_value(x.start);
  boost::hash_combine(result, x.end);
  return result;
}
}  // namespace majiq
namespace std {
// specialize std::hash for intervals
template <> struct hash<majiq::OpenInterval> {
  std::size_t operator()(const majiq::OpenInterval& x) const noexcept {
    return majiq::hash_value(x);
  }
};
template <> struct hash<majiq::ClosedInterval> {
  std::size_t operator()(const majiq::ClosedInterval& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std
#endif  // MAJIQ_INTERVAL_HPP
