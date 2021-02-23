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
#include <algorithm>
#include <numeric>
#include <utility>

#include "SJJunctions.hpp"
#include "Exons.hpp"
#include "GeneIntrons.hpp"
#include "ContigIntrons.hpp"
#include "MajiqConstants.hpp"


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
  const junction_pos_t num_positions_;

  void check_valid() const;

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
        num_positions_{num_positions} { check_valid(); }
  SJJunctionsPositions(
      std::shared_ptr<SJJunctions>&& junctions,
      std::vector<PositionReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t num_positions)
      : junctions_{junctions},
        reads_{reads},
        offsets_{offsets},
        num_positions_{num_positions} { check_valid(); }
  SJJunctionsPositions(const SJJunctionsPositions& x) = default;
  SJJunctionsPositions(SJJunctionsPositions&& x) = default;
  SJJunctionsPositions& operator=(const SJJunctionsPositions& x) = delete;
  SJJunctionsPositions& operator=(SJJunctionsPositions&& x) = delete;
};

namespace detail {
/**
 * Define how raw positions are coarsely binned to match junctions
 */
struct IntronCoarseBins {
  // how many positions in smallest bin
  const junction_pos_t min_pos_per_bin;
  // how many bins with exactly 1 more than min
  const junction_pos_t bins_with_extra;

  static junction_pos_t num_raw_positions(
      junction_pos_t intron_length, junction_pos_t num_bins) {
    return intron_length + num_bins;
  }

  IntronCoarseBins(junction_pos_t num_raw_positions, junction_pos_t num_bins)
      : min_pos_per_bin{num_raw_positions / num_bins},
        bins_with_extra{num_raw_positions % num_bins} { }

  junction_pos_t resulting_bin(junction_pos_t raw_read_position) const {
    // bins with extra position go in front, find position that divides that
    const junction_pos_t first_pos_without_extra
      = bins_with_extra * (1 + min_pos_per_bin);
    // determine where to put position relative to that
    return (
        raw_read_position < first_pos_without_extra
        ? raw_read_position / (1 + min_pos_per_bin)
        : (bins_with_extra
          + ((raw_read_position - first_pos_without_extra)
            / min_pos_per_bin)));
  }
  junction_pos_t bin_num_positions(junction_pos_t bin_idx) const {
    return bin_idx < bins_with_extra ? 1 + min_pos_per_bin : min_pos_per_bin;
  }
};
}  // namespace detail

class SJIntronsBins {
 private:
  const ContigIntrons introns_;
  const std::vector<PositionReads> reads_;
  const std::vector<size_t> offsets_;
  const junction_pos_t num_bins_;

  void check_valid() const;  // raise error if defined inconsistently

 public:
  static SJIntronsBins FromBam(
      const char* infile, const junction_pos_t num_bins, const Exons& exons,
      const GeneIntrons& gene_introns, ExperimentStrandness exp_strandness,
      int nthreads);

  size_t num_introns() const noexcept { return introns_.size(); }
  size_t size() const noexcept { return reads_.size(); }
  junction_pos_t num_bins() const noexcept { return num_bins_; }
  const ContigIntrons& introns() { return introns_; }
  const std::vector<PositionReads>& reads() { return reads_; }
  const std::vector<size_t>& offsets() { return offsets_; }

  // summaries per intron
  junction_pos_t numbins_nonzero(size_t intron_idx) const {
    return offsets_[1 + intron_idx] - offsets_[intron_idx];
  }
  junction_pos_t numbins_mincov(size_t intron_idx, junction_ct_t mincov) const {
    if (mincov == 0) {
      return numbins_nonzero(intron_idx);
    } else {
      return std::count_if(
          reads_.begin() + offsets_[intron_idx],
          reads_.begin() + offsets_[1 + intron_idx],
          [&mincov](const PositionReads& x) { return x.reads >= mincov; });
    }
  }
  junction_ct_t numreads(size_t intron_idx, junction_pos_t num_stacks) const {
    size_t nonstack_end
      = num_stacks == 0 || num_stacks < numbins_nonzero(intron_idx)
      ? offsets_[1 + intron_idx] - num_stacks : offsets_[intron_idx];
    return std::accumulate(
        reads_.begin() + offsets_[intron_idx],
        reads_.begin() + nonstack_end,
        junction_ct_t{},
        [](junction_ct_t s, const PositionReads& x) { return s + x.reads; });
  }
  real_t scaled_numreads(size_t intron_idx, junction_pos_t num_stacks) const {
    using detail::IntronCoarseBins;
    auto intron_length = introns_[intron_idx].coordinates.length();
    return static_cast<real_t>(numreads(intron_idx, num_stacks) * intron_length)
      / IntronCoarseBins::num_raw_positions(intron_length, num_bins_);
  }
  const PositionReads& reads_elem(
      size_t intron_idx, junction_pos_t nonzero_idx) const {
    // NOTE: assumes nonzero_idx < numbins_nonzero(intron_idx)
    return reads_[offsets_[intron_idx] + nonzero_idx];
  }
  // scaled NOT to junctions but to fractional positions per bin (average)
  real_t scaled_reads_elem(
      size_t intron_idx, junction_pos_t nonzero_idx) const {
    using detail::IntronCoarseBins;
    // get average number of positions per bin
    auto total_positions = IntronCoarseBins::num_raw_positions(
        introns_[intron_idx].coordinates.length(), num_bins_);
    real_t avg_num_positions = static_cast<real_t>(total_positions) / num_bins_;
    // compare to actual element
    const auto& raw_elem = reads_elem(intron_idx, nonzero_idx);
    auto bin_num_positions = (
        IntronCoarseBins(total_positions, num_bins_)
        .bin_num_positions(raw_elem.pos));
    // rescale
    return raw_elem.reads * avg_num_positions / bin_num_positions;
  }

  SJIntronsBins(
      ContigIntrons&& introns,
      std::vector<PositionReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t num_bins)
      : introns_{std::move(introns)},
        reads_{std::move(reads)},
        offsets_{std::move(offsets)},
        num_bins_{num_bins} { check_valid(); }
  SJIntronsBins(const SJIntronsBins&) = default;
  SJIntronsBins(SJIntronsBins&&) = default;
  SJIntronsBins& operator=(const SJIntronsBins&) = delete;
  SJIntronsBins& operator=(SJIntronsBins&&) = delete;
};

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONSPOSITIONS_HPP
