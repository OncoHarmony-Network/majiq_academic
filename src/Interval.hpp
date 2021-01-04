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

#include "MajiqTypes.hpp"

namespace majiq {
namespace detail {

struct Interval {
 public:
  position_t start;
  position_t end;

  // constructors
  Interval(position_t a, position_t b) : start{a}, end{b} { }
  Interval(const Interval& x) : Interval(x.start, x.end) { }

  // full interval or otherwise?
  bool has_start() const { return start >= 0; }
  bool has_end() const { return end >= 0; }
  bool is_full_interval() const { return has_start() && has_end(); }
  bool is_half_interval() const { return has_start() != has_end(); }
  bool is_invalid() const { return !(has_start() || has_end()); }

  // valid positions (if any)
  const position_t& first_pos() const { return has_start() ? start : end; }
  const position_t& last_pos() const { return has_end() ? end : start; }

  // ordering
  bool operator<(const Interval& rhs) const {
    return std::tie(first_pos(), last_pos())
      < std::tie(rhs.first_pos(), rhs.last_pos());
  }
  // equality
  bool operator==(const Interval& rhs) const {
    return start == rhs.start && end == rhs.end;
  }
};

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
      throw std::logic_error("ClosedInterval must have nonnegative length");
    }
  }
  ClosedInterval(const ClosedInterval& x) : ClosedInterval(x.start, x.end) { }
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
      throw std::logic_error("OpenInterval must have nonnegative length");
    }
  }
  OpenInterval(const OpenInterval& x) : OpenInterval(x.start, x.end) { }
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
std::ostream& operator<<(std::ostream& os, const ClosedInterval& x) {
  os << "[" << x.start << ", " << x.end << "]";
  return os;
}
std::ostream& operator<<(std::ostream& os, const OpenInterval& x) {
  os << "(" << x.start << ", " << x.end << ")";
  return os;
}


}  // namespace majiq
#endif  // MAJIQ_INTERVAL_HPP
