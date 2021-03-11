/**
 * MaijqTypes.hpp
 *
 * base types for MAJIQ
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_TYPES_HPP
#define MAJIQ_TYPES_HPP

#include <cstdint>
#include <string>
#include <iostream>
#include <functional>
#include <random>

namespace majiq {

struct EmptyDataT {};

using real_t = float;  // for numerics
using position_t = int64_t;  // for coordinates, difference in coordinates
using junction_pos_t = int32_t;  // position for junctions
using junction_ct_t = int32_t;  // counts of reads with splits
using intron_ct_t = real_t;  // count for introns

using seqid_t = std::string;  // type for contig ids
using geneid_t = std::string;
using genename_t = std::string;

using rng_t = std::mt19937;

enum class GeneStrandness : unsigned char {
  FORWARD = '+',
  REVERSE = '-',
  AMBIGUOUS = '.',
};
inline std::ostream& operator<<(
    std::ostream& os, const GeneStrandness& x) noexcept {
  os << static_cast<char>(x);
  return os;
}

enum class ExperimentStrandness : unsigned char {
  FORWARD = 'F',  // read 1 is in forward direction, Salmon ISF
  REVERSE = 'R',  // read 2 is in forward direction, Salmon ISR
  NONE = 'N',  // could be either way, Salmon IU
};

// did a junction pass different thresholds?
enum class JunctionPassedStatus : unsigned char {
  NOT_PASSED,
  ANNOTATED_PASSED,
  DENOVO_PASSED
};

enum class EventType : unsigned char {
  SRC_EVENT = 's',
  DST_EVENT = 't'
};
inline EventType OtherEventType(EventType x) {
  switch (x) {
    case EventType::SRC_EVENT:
      return EventType::DST_EVENT;
    case EventType::DST_EVENT:
      return EventType::SRC_EVENT;
    default:
      return x;
  }
}
inline std::ostream& operator<<(
    std::ostream& os, const EventType& x) noexcept {
  os << static_cast<char>(x);
  return os;
}

struct Event {
  size_t ref_exon_idx_;
  EventType type_;
  Event() = default;
  Event(size_t ref_exon_idx, EventType type)
      : ref_exon_idx_{ref_exon_idx}, type_{type} { }
};
inline bool operator<(const Event& x, const Event& y) noexcept {
  return std::tie(x.ref_exon_idx_, x.type_)
    < std::tie(y.ref_exon_idx_, y.type_);
}
inline bool operator>(const Event& x, const Event& y) noexcept {
  return y < x;
}
inline bool operator>=(const Event& x, const Event& y) noexcept {
  return !(x < y);
}
inline bool operator<=(const Event& x, const Event& y) noexcept {
  return !(y < x);
}
inline bool operator==(const Event& x, const Event& y) noexcept {
  return std::tie(x.ref_exon_idx_, x.type_)
    == std::tie(y.ref_exon_idx_, y.type_);
}
inline bool operator!=(const Event& x, const Event& y) noexcept {
  return !(x == y);
}
}  // namespace majiq


#endif  // MAJIQ_TYPES_HPP
