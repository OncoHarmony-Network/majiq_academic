/**
 * SJJunctionsPositions.hpp
 *
 * Underlying raw per-position coverage for SJJunctions
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_SJJUNCTIONSPOSITIONS_HPP
#define MAJIQ_SJJUNCTIONSPOSITIONS_HPP

#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <stdexcept>

#include "SJJunctions.hpp"


namespace majiq {
// identity/read count for nonzero junction positions
struct PositionReads {
  junction_pos_t pos;
  junction_ct_t reads;
};
inline bool operator<(const PositionReads& x, const PositionReads& y) noexcept {
  // sorted by reads first
  return std::tie(x.reads, x.pos) < std::tie(y.reads, y.pos);
}

class SJJunctionsPositions {
 private:
  std::shared_ptr<SJJunctions> junctions_;
  std::vector<PositionReads> reads_;
  std::vector<size_t> offsets_;
  junction_pos_t num_positions_;

  // check that input values together are valid
  inline static bool is_valid(
      const std::shared_ptr<SJJunctions>& junctions,
      const std::vector<PositionReads>& reads,
      const std::vector<size_t>& offsets,
      junction_pos_t num_positions);

 public:
  inline size_t num_junctions() const noexcept { return junctions_->size(); }
  inline size_t size() const noexcept { return reads_.size(); }
  inline junction_pos_t num_positions() const { return num_positions_; }
  const std::shared_ptr<SJJunctions>& junctions() { return junctions_; }
  const std::vector<PositionReads>& reads() { return reads_; }
  const std::vector<size_t>& offsets() { return offsets_; }

  SJJunctionsPositions()
      : junctions_{std::make_shared<SJJunctions>()},
        reads_{},
        offsets_{0},
        num_positions_{0} { }
  SJJunctionsPositions(
      const std::shared_ptr<SJJunctions>& junctions,
      const std::vector<PositionReads>& reads,
      const std::vector<size_t>& offsets,
      junction_pos_t num_positions)
      : junctions_{junctions},
        reads_{reads},
        offsets_{offsets},
        num_positions_{num_positions} {
    if (!is_valid(junctions_, reads_, offsets_, num_positions_)) {
      throw std::invalid_argument(
          "SJJunctionsPositions given invalid arguments");
    }
  }
  SJJunctionsPositions(
      std::shared_ptr<SJJunctions>&& junctions,
      std::vector<PositionReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t num_positions)
      : junctions_{junctions},
        reads_{reads},
        offsets_{offsets},
        num_positions_{num_positions} {
    if (!is_valid(junctions_, reads_, offsets_, num_positions_)) {
      throw std::invalid_argument(
          "SJJunctionsPositions given invalid arguments");
    }
  }
  SJJunctionsPositions(const SJJunctionsPositions& x) = default;
  SJJunctionsPositions(SJJunctionsPositions&& x) = default;
  SJJunctionsPositions& operator=(const SJJunctionsPositions& x) = default;
  SJJunctionsPositions& operator=(SJJunctionsPositions&& x) = default;
};

inline bool SJJunctionsPositions::is_valid(
    const std::shared_ptr<SJJunctions>& junctions,
    const std::vector<PositionReads>& reads,
    const std::vector<size_t>& offsets,
    junction_pos_t num_positions) {
  if (
      // must have junctions
      junctions == nullptr
      // offsets[i], offsets[i + 1] correspond to junctions[i] --> match sizes
      || junctions->size() != offsets.size() - 1
      // last offset must correspond to size of reads
      || offsets.back() != reads.size()) {
    return false;
  }
  // for each junction
  for (size_t jidx = 0; jidx < junctions->size(); ++jidx) {
    // offsets into junction positions
    auto jp_start = reads.begin() + offsets[jidx];
    auto jp_end = reads.begin() + offsets[jidx + 1];
    if (
        // offsets must be non-decreasing
        jp_end < jp_start
        // difference in offsets correspond to numpos
        || jp_end - jp_start != (*junctions)[jidx].data.numpos
        // can't have numpos greater than total possible
        || (*junctions)[jidx].data.numpos > num_positions
        // all junction positions must be in [0, num_positions)
        || !std::all_of(jp_start, jp_end,
          [num_positions](const PositionReads& x) {
          return x.pos < num_positions;
          })
        // sum of position reads equals junction numreads
        || std::accumulate(jp_start, jp_end, junction_ct_t{},
          [](junction_ct_t s, const PositionReads& x) { return s + x.reads; })
          != (*junctions)[jidx].data.numreads
        // PositionReads in sorted order (by number of reads)
        || !std::is_sorted(jp_start, jp_end)) {
      return false;
    }
  }
  return true;
}

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONSPOSITIONS_HPP
