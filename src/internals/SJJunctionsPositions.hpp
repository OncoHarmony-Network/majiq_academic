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
// identity/read count for nonzero intron/junction bins
struct BinReads {
  junction_pos_t bin_idx;
  junction_ct_t bin_reads;
};
inline bool operator<(const BinReads& x, const BinReads& y) noexcept {
  return std::tie(x.bin_reads, x.bin_idx) < std::tie(y.bin_reads, y.bin_idx);
}

namespace detail {
template <typename RegionsT>
class SJRegionBinReads {
 protected:
  const std::shared_ptr<RegionsT> regions_;
  const std::vector<BinReads> reads_;
  const std::vector<size_t> offsets_;
  const junction_pos_t total_bins_;

 public:
  size_t num_regions() const noexcept { return regions_->size(); }
  size_t size() const noexcept { return reads_.size(); }
  junction_pos_t total_bins() const noexcept { return total_bins_; }
  const std::shared_ptr<RegionsT>& regions() { return regions_; }
  const std::vector<BinReads>& reads() { return reads_; }
  const std::vector<size_t>& offsets() { return offsets_; }

  // per region summaries
  junction_pos_t numbins_nonzero(size_t i) const {
    return offsets_[1 + i] - offsets_[i];
  }
  junction_pos_t numbins_minreads(size_t i, junction_ct_t minreads) const {
    if (minreads == 0) {
      return numbins_nonzero(i);
    } else {
      return std::count_if(
          reads_.begin() + offsets_[i], reads_.begin() + offsets_[1 + i],
          [minreads](const BinReads& x) { return x.bin_reads >= minreads; });
    }
  }
  junction_ct_t numreads(size_t i, junction_pos_t num_stacks) const {
    const size_t nonstack_end
      = num_stacks == 0 ? offsets_[1 + i] : (
          num_stacks < numbins_nonzero(i)
          ? offsets_[1 + i] - num_stacks : offsets_[i]);
    return std::accumulate(
        reads_.begin() + offsets_[i], reads_.begin() + nonstack_end,
        junction_ct_t{},
        [](junction_ct_t s, const BinReads& x) { return s + x.bin_reads; });
  }
  const BinReads& reads_elem(size_t i, junction_pos_t nonzero_idx) const {
    return reads_[offsets_[i] + nonzero_idx];
  }

 private:
  void check_valid() const {
    if (regions_ == nullptr) {
      throw std::runtime_error("SJRegionBinReads object has null regions");
    } else if (regions_->size() + 1 != offsets_.size()) {
      throw std::runtime_error(
          "SJRegionBinReads offsets do not correspond to regions");
    } else if (offsets_.back() != reads_.size()) {
      throw std::runtime_error(
          "SJRegionBinReads offsets do not correspond to reads per bin");
    } else {
      for (size_t i = 0; i < regions_->size(); ++i) {
        if (offsets_[i] > offsets_[1 + i]) {
          throw std::runtime_error(
              "SJRegionBinReads has nonmonotonically increasing offsets");
        } else if (numbins_nonzero(i) > total_bins_) {
          throw std::runtime_error(
              "SJRegionBinReads has region with more bins than possible");
        } else if (
            !std::all_of(
              reads_.begin() + offsets_[i], reads_.begin() + offsets_[1 + i],
              [total_bins = total_bins_](const BinReads& x) {
              return x.bin_idx < total_bins; })) {
          throw std::runtime_error(
              "SJRegionBinReads has reported bin_idx outside valid range");
        } else if (!std::is_sorted(
              reads_.begin() + offsets_[i], reads_.begin() + offsets_[1 + i])) {
          throw std::runtime_error(
              "SJRegionBinReads is not appropriately sorted per region");
        }
      }  // loop over regions
    }
  }

 public:
  SJRegionBinReads(
      const std::shared_ptr<RegionsT>& regions,
      std::vector<BinReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t total_bins)
      : regions_{regions},
        reads_{std::move(reads)},
        offsets_{std::move(offsets)},
        total_bins_{total_bins} { }
  SJRegionBinReads(const SJRegionBinReads&) = default;
  SJRegionBinReads(SJRegionBinReads&&) = default;
  SJRegionBinReads& operator=(const SJRegionBinReads&) = delete;
  SJRegionBinReads& operator=(SJRegionBinReads&&) = delete;
};  // class SJRegionPositions

}  // namespace detail

class SJJunctionsPositions : public detail::SJRegionBinReads<SJJunctions> {
  using BaseT = detail::SJRegionBinReads<SJJunctions>;

 private:
  // already have check_valid from parent class, but check SJJunctions data
  void check_junctions() const {
    for (size_t i = 0; i < regions_->size(); ++i) {
      const SJJunction& curj = (*regions_)[i];
      if (curj.numpos() != numbins_nonzero(i)) {
        throw std::runtime_error("SJJunctions/Positions mismatch on numpos");
      } else if (curj.numreads() != numreads(i, 0)) {
        throw std::runtime_error("SJJunctions/Positions mismatch on numreads");
      }
    }
  }

 public:
  static SJJunctionsPositions FromBam(
      const char* infile, ExperimentStrandness exp_strandness, int nthreads);

  SJJunctionsPositions(
      const std::shared_ptr<SJJunctions>& junctions,
      std::vector<BinReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t num_positions)
      : BaseT{junctions, std::move(reads), std::move(offsets), num_positions} {
    check_junctions();
  }
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

class SJIntronsBins : public detail::SJRegionBinReads<ContigIntrons> {
  using BaseT = detail::SJRegionBinReads<ContigIntrons>;

 public:
  static SJIntronsBins FromBam(
      const char* infile, const junction_pos_t num_bins, const Exons& exons,
      const GeneIntrons& gene_introns, ExperimentStrandness exp_strandness,
      int nthreads);

  real_t scaled_numreads(size_t i, junction_pos_t num_stacks) const {
    using detail::IntronCoarseBins;
    const auto intron_length = (*regions_)[i].coordinates.length();
    return static_cast<real_t>(numreads(i, num_stacks) * intron_length)
      / IntronCoarseBins::num_raw_positions(intron_length, total_bins_);
  }
  real_t scaled_bin_reads_elem(size_t i, junction_pos_t nonzero_idx) const {
    using detail::IntronCoarseBins;
    // get average number of positions per bin
    const auto total_positions = IntronCoarseBins::num_raw_positions(
        (*regions_)[i].coordinates.length(), total_bins_);
    const real_t avg_num_positions
      = static_cast<real_t>(total_positions) / total_bins_;
    // compare to actual element
    const auto& raw_elem = reads_elem(i, nonzero_idx);
    const auto bin_num_positions = (
        IntronCoarseBins(total_positions, total_bins_)
        .bin_num_positions(raw_elem.bin_idx));
    // rescale
    return raw_elem.bin_reads * avg_num_positions / bin_num_positions;
  }

  SJIntronsBins(
      const std::shared_ptr<ContigIntrons>& introns,
      std::vector<BinReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t total_bins)
      : BaseT{introns, std::move(reads), std::move(offsets), total_bins} { }
  SJIntronsBins(const SJIntronsBins& x) = default;
  SJIntronsBins(SJIntronsBins&& x) = default;
  SJIntronsBins& operator=(const SJIntronsBins& x) = delete;
  SJIntronsBins& operator=(SJIntronsBins&& x) = delete;
};

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONSPOSITIONS_HPP
