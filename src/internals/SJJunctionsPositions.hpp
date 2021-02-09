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

  // used by FromBam
  static constexpr uint32_t USE_MIN_OVERHANG = 8;

  // check that input values together are valid
  static bool is_valid(
      const std::shared_ptr<SJJunctions>& junctions,
      const std::vector<PositionReads>& reads,
      const std::vector<size_t>& offsets,
      junction_pos_t num_positions);

 public:
  static SJJunctionsPositions FromBam(
      const char* infile, ExperimentStrandness exp_strandness, int nthreads);

  size_t num_junctions() const noexcept { return junctions_->size(); }
  size_t size() const noexcept { return reads_.size(); }
  junction_pos_t num_positions() const { return num_positions_; }
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

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONSPOSITIONS_HPP
