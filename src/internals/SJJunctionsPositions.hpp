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
#include <random>

#include <boost/math/distributions/poisson.hpp>

#include "MajiqTypes.hpp"
#include "SJJunctions.hpp"
#include "Exons.hpp"
#include "GeneIntrons.hpp"
#include "SJIntrons.hpp"
#include "MajiqConstants.hpp"
#include "ExperimentThresholds.hpp"


namespace majiq {
// identity/read count for nonzero intron/junction bins
template <typename CountT>
struct BinReads {
  junction_pos_t bin_idx;
  CountT bin_reads;
};
template <typename T>
inline bool operator<(const BinReads<T>& x, const BinReads<T>& y) noexcept {
  return std::tie(x.bin_reads, x.bin_idx) < std::tie(y.bin_reads, y.bin_idx);
}

namespace detail {
template <typename RegionsT, typename CountT>
class SJRegionBinReads {
 protected:
  const std::shared_ptr<RegionsT> regions_;
  const std::vector<BinReads<CountT>> reads_;
  const std::vector<size_t> offsets_;
  const junction_pos_t total_bins_;

 public:
  size_t num_regions() const noexcept { return regions_->size(); }
  size_t size() const noexcept { return reads_.size(); }
  junction_pos_t total_bins() const noexcept { return total_bins_; }
  const std::shared_ptr<RegionsT>& regions() const { return regions_; }
  const std::vector<BinReads<CountT>>& reads() { return reads_; }
  const std::vector<size_t>& offsets() { return offsets_; }

  // per region summaries
  junction_pos_t numbins_nonzero(size_t i) const {
    return offsets_[1 + i] - offsets_[i];
  }

  typename std::vector<BinReads<CountT>>::const_iterator
  begin_region(size_t i) const { return reads_.cbegin() + offsets_[i]; }
  typename std::vector<BinReads<CountT>>::const_iterator
  end_region(size_t i) const { return begin_region(1 + i); }
  typename std::vector<BinReads<CountT>>::const_iterator
  end_region(size_t i, junction_pos_t num_stacks) const {
    return num_stacks == 0
      ? end_region(i)
      : (num_stacks < numbins_nonzero(i)
          ? end_region(i) - num_stacks : begin_region(i));
  }

  junction_pos_t numbins_minreads(size_t i, CountT minreads) const {
    if (minreads == 0) {
      return numbins_nonzero(i);
    } else {
      auto end = end_region(i);
      return end - std::find_if(begin_region(i), end,
          [minreads](const BinReads<CountT>& x) {
          return x.bin_reads >= minreads;
          });
    }
  }
  CountT numreads(size_t i, junction_pos_t num_stacks) const {
    return std::accumulate(
        begin_region(i), end_region(i, num_stacks),
        CountT{},
        [](CountT s, const BinReads<CountT>& x) { return s + x.bin_reads; });
  }
  const BinReads<CountT>& reads_elem(
      size_t i, junction_pos_t nonzero_idx) const {
    return *(begin_region(i) + nonzero_idx);
  }

  template <typename OutputIt>
  void bootstrap_coverage(
      size_t i, junction_pos_t num_stacks,
      rng_t& generator,
      OutputIt first, OutputIt last) const {
    using OutputT = remove_cvref_t<decltype(*first)>;
    const auto pos_begin = begin_region(i);
    const auto pos_end = end_region(i, num_stacks);
    if (pos_begin == pos_end) {
      std::fill(first, last, OutputT{0});
      return;
    }
    // do we resort to parametric sampling? We know if all values are same, yes
    bool parametric_bootstrap
      = pos_begin->bin_reads == (pos_end - 1)->bin_reads;
    const junction_pos_t num_nonzero = pos_end - pos_begin;
    OutputT bootstrap_mean;
    // otherwise, we have to check empirical mean/variance of nonzero positions
    if (parametric_bootstrap) {
      bootstrap_mean = static_cast<OutputT>(pos_begin->bin_reads);
    } else {
      bootstrap_mean = std::accumulate(
          pos_begin, pos_end, OutputT{0},
          [](OutputT s, const BinReads<CountT>& x) {
          return s + static_cast<OutputT>(x.bin_reads); })
        / num_nonzero;
      OutputT bootstrap_variance = std::accumulate(
          pos_begin, pos_end, OutputT{0},
          [bootstrap_mean](OutputT s, const BinReads<CountT>& x) {
          OutputT residual = static_cast<OutputT>(x.bin_reads - bootstrap_mean);
          return s + residual * residual; })
        / num_nonzero;  // population variance, no correction or sample variance
      parametric_bootstrap = bootstrap_variance < bootstrap_mean;
    }
    // perform the bootstrapping
    if (parametric_bootstrap) {
      // bootstrap distribution is Poisson in parametric case
      std::poisson_distribution<int> dist{bootstrap_mean * num_nonzero};
      std::generate(first, last,
          [&dist, &generator]() -> OutputT {
          return static_cast<OutputT>(dist(generator)); });
    } else {
      // sample one less than nonzero positions with replacement, rescale
      const junction_pos_t num_samples = num_nonzero - 1;
      std::uniform_int_distribution<junction_pos_t> dist{0, num_nonzero - 1};
      auto random_bin_reads = [&dist, &generator, &pos_begin]() {
        return pos_begin[dist(generator)].bin_reads; };
      auto sample = [&num_nonzero, &num_samples, &random_bin_reads]() {
        OutputT result = 0;
        for (junction_pos_t i = 0; i < num_samples; ++i) {
          result += random_bin_reads();
        }
        return static_cast<OutputT>(result * num_nonzero) / num_samples; };
      std::generate(first, last, sample);
    }
  }

  junction_pos_t numstacks(size_t i, real_t pvalue_threshold) const {
    using dist_t = boost::math::poisson_distribution<real_t>;
    using boost::math::cdf;
    using boost::math::complement;
    const junction_pos_t numbins = numbins_nonzero(i);
    if (numbins == 0) {
      return 0;
    } else if (numbins == 1) {
      // test against Poisson rate 0.5
      constexpr real_t SINGLETON_RATE_H0 = 0.5;
      const dist_t stack_dist(SINGLETON_RATE_H0);
      return (
          cdf(complement(stack_dist, reads_elem(i, 0).bin_reads))
          < pvalue_threshold)
        ? 1 : 0;
    } else {  // numbins > 1
      // calculate by comparing to leave-one-out mean
      // TODO(jaicher): consider replacing with median (or variants of median)
      CountT sum = numreads(i, 0);
      auto it = end_region(i);  // iterate backwards from here
      do {
        --it;
        real_t loo_mean
          = static_cast<real_t>(sum - it->bin_reads) / (numbins - 1);
        const dist_t stack_dist(loo_mean);
        if (cdf(complement(stack_dist, it->bin_reads)) >= pvalue_threshold) {
          // it is not a stack, and so are any of the bins before it
          break;
        }
      } while (it != begin_region(i));
      return (end_region(i) - it) - 1;
    }
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
              [total_bins = total_bins_](const BinReads<CountT>& x) {
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
      std::vector<BinReads<CountT>>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t total_bins)
      : regions_{regions},
        reads_{std::move(reads)},
        offsets_{std::move(offsets)},
        total_bins_{total_bins} {
    if (regions_ == nullptr) {
      throw std::invalid_argument("SJRegionBinReads requires non-null regions");
    }
  }
};  // class SJRegionPositions

}  // namespace detail

class SJJunctionsPositions
    : public detail::SJRegionBinReads<SJJunctions, junction_ct_t> {
  using BaseT = detail::SJRegionBinReads<SJJunctions, junction_ct_t>;

 public:
  junction_pos_t weight(size_t i) const { return total_bins(); }
  static SJJunctionsPositions FromBam(
      const char* infile, ExperimentStrandness exp_strandness, int nthreads);

  JunctionPassedStatus passed(
      size_t i, const ExperimentThresholds& thresholds) const {
    junction_pos_t numpos = numbins_nonzero(i);
    if (numpos < thresholds.minpos_) {
      return JunctionPassedStatus::NOT_PASSED;
    } else {
      constexpr junction_pos_t NOSTACKS = 0;  // raw reads
      junction_ct_t numreads = this->numreads(i, NOSTACKS);
      if (numreads < thresholds.minreads_) {
        return JunctionPassedStatus::NOT_PASSED;
      } else if (numreads < thresholds.mindenovo_) {
        return JunctionPassedStatus::ANNOTATED_PASSED;
      } else {
        return JunctionPassedStatus::DENOVO_PASSED;
      }
    }
  }

  SJJunctionsPositions(
      const std::shared_ptr<SJJunctions>& junctions,
      std::vector<BinReads<junction_ct_t>>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t num_positions)
      : BaseT{junctions, std::move(reads), std::move(offsets), num_positions} {
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

class SJIntronsBins
    : public detail::SJRegionBinReads<SJIntrons, intron_ct_t> {
  using BaseT = detail::SJRegionBinReads<SJIntrons, intron_ct_t>;

 public:
  static SJIntronsBins FromBam(
      const char* infile, const junction_pos_t num_bins, const Exons& exons,
      const GeneIntrons& gene_introns, ExperimentStrandness exp_strandness,
      int nthreads);

  junction_pos_t num_raw_positions(size_t i) const {
    using detail::IntronCoarseBins;
    const auto intron_length = (*regions_)[i].coordinates.length();
    return IntronCoarseBins::num_raw_positions(intron_length, total_bins());
  }
  junction_pos_t weight(size_t i) const { return num_raw_positions(i); }

  real_t scaled_numreads(size_t i, junction_pos_t num_stacks) const {
    const auto intron_length = (*regions_)[i].coordinates.length();
    return static_cast<real_t>(numreads(i, num_stacks) * intron_length)
      / num_raw_positions(i);
  }
  bool passed(size_t i, const IntronThresholdsGenerator& it_gen,
      junction_pos_t num_stacks) const {
    // get thresholds for this intron
    IntronThresholds thresholds = it_gen((*regions_)[i].coordinates.length());
    return scaled_numreads(i, num_stacks) >= thresholds.minreads_
      && (numbins_minreads(i, thresholds.mincov_)
          >= num_stacks + thresholds.minbins_);
  }

  SJIntronsBins(
      const std::shared_ptr<SJIntrons>& introns,
      std::vector<BinReads<intron_ct_t>>&& reads,
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
