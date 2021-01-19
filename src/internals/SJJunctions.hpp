/**
 * SJJunctions.hpp
 *
 * Contig regions with different data to represent junctions in different
 * settings
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_SJJUNCTIONS_HPP
#define MAJIQ_SJJUNCTIONS_HPP

#include <vector>
#include <map>
#include <utility>
#include <memory>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <string>

#include "MajiqTypes.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"

namespace majiq {
// total reads and nonzero positions for a junction
struct ExperimentCounts {
  junction_ct_t numreads;
  junction_pos_t numpos;
};
// group passed (build, denovo)
struct GroupPassed {
  bool build_passed;
  bool denovo_passed;
};


using ContigJunction = detail::ContigRegion<OpenInterval, detail::EmptyDataT>;
using SJJunction = detail::ContigRegion<OpenInterval, ExperimentCounts>;
using GroupJunction = detail::ContigRegion<OpenInterval, GroupPassed>;

using SJJunctions = detail::ContigRegions<SJJunction>;
using GroupJunctions = detail::ContigRegions<GroupJunction>;


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

  static bool is_valid(
      const std::shared_ptr<SJJunctions>& junctions,
      const std::vector<PositionReads>& reads,
      const std::vector<size_t>& offsets,
      junction_pos_t num_positions) {
    if (junctions == nullptr
        || offsets.empty()
        || junctions->size() != offsets.size() - 1
        || offsets.back() != reads.size()) {
      return false;
    }
    for (size_t jidx = 0; jidx < junctions->size(); ++jidx) {
      auto jp_start = reads.begin() + offsets[jidx];
      auto jp_end = reads.begin() + offsets[jidx + 1];
      if (jp_end < jp_start
          || jp_end - jp_start != (*junctions)[jidx].data.numpos
          || (*junctions)[jidx].data.numpos > num_positions
          || !std::all_of(jp_start, jp_end,
            [num_positions](const PositionReads& x) {
            return x.pos < num_positions;
            })
          || std::accumulate(jp_start, jp_end, junction_ct_t{},
            [](junction_ct_t s, const PositionReads& x) { return s + x.reads; })
            != (*junctions)[jidx].data.numreads
          || !std::is_sorted(jp_start, jp_end)) {
        return false;
      }
    }
    return true;
  }

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
      junction_pos_t num_positions) {
    if (!is_valid(junctions, reads, offsets, num_positions)) {
      throw std::invalid_argument(
          "SJJunctionsPositions given invalid arguments");
    }
    junctions_ = junctions;
    reads_ = reads;
    offsets_ = offsets;
    num_positions_ = num_positions;
  }
  SJJunctionsPositions(
      std::shared_ptr<SJJunctions>&& junctions,
      std::vector<PositionReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t num_positions) {
    if (!is_valid(junctions, reads, offsets, num_positions)) {
      throw std::invalid_argument(
          "SJJunctionsPositions given invalid arguments");
    }
    junctions_ = junctions;
    reads_ = reads;
    offsets_ = offsets;
    num_positions_ = num_positions;
  }
  SJJunctionsPositions(const SJJunctionsPositions& x) = default;
  SJJunctionsPositions(SJJunctionsPositions&& x) = default;
  SJJunctionsPositions& operator=(const SJJunctionsPositions& x) = default;
  SJJunctionsPositions& operator=(SJJunctionsPositions&& x) = default;
};


// how do we get from SJJunctions to GroupJunctions?
class GroupJunctionsGenerator {
 private:
  const junction_ct_t minreads_;
  const junction_ct_t mindenovo_;  // REQUIRED >= minreads
  const junction_pos_t minpos_;

  std::shared_ptr<Contigs> contigs_;

  // number of experiments passing filters
  struct GroupCounts {
    uint32_t numpassed;  // passing build filters
    uint32_t numdenovo;  // passing denovo filters
  };
  size_t num_experiments_;
  std::vector<  // over contigs
    std::map<std::pair<OpenInterval, GeneStrandness>, GroupCounts>> counts_;

 public:
  GroupJunctionsGenerator(junction_ct_t minreads, junction_ct_t mindenovo,
      junction_pos_t minpos, std::shared_ptr<Contigs> contigs)
      : minreads_{minreads}, mindenovo_{mindenovo}, minpos_{minpos},
        contigs_{contigs == nullptr ? std::make_shared<Contigs>() : contigs} {
    if (minreads > mindenovo) {
      throw std::invalid_argument("mindenovo must be at least minreads");
    }
  }
  GroupJunctionsGenerator(junction_ct_t minreads, junction_ct_t mindenovo,
      junction_pos_t minpos)
      : GroupJunctionsGenerator{minreads, mindenovo, minpos, nullptr} {
  }
  GroupJunctionsGenerator(const GroupJunctionsGenerator& x) = default;
  GroupJunctionsGenerator(GroupJunctionsGenerator&& x) = default;
  GroupJunctionsGenerator& operator=(
      const GroupJunctionsGenerator& x) = delete;
  GroupJunctionsGenerator& operator=(GroupJunctionsGenerator&& x) = delete;

  /**
   * Create group junctions with data aggregated from input experiments,
   * filtering unpassed junctions using min_experiment
   */
  GroupJunctions PassedJunctions(size_t min_experiments) {
    std::vector<GroupJunction> passed;
    for (size_t contig_idx = 0; contig_idx < contigs_->size(); ++contig_idx) {
      for (auto&& [coord_strand, ct] : counts_[contig_idx]) {
        if (ct.numpassed >= min_experiments) {
          passed.emplace_back(
              (*contigs_)[contig_idx], coord_strand.first, coord_strand.second,
              GroupPassed{ct.numpassed >= min_experiments,
                          ct.numdenovo >= min_experiments});
        }
      }
    }
    // create GroupJunctions
    return GroupJunctions{passed};
  }
  /**
   * Add junctions from experiment to group towards group filters
   */
  void AddExperiment(const SJJunctions& experiment) {
    for (const SJJunction& junction : experiment) {
      // need min positions to pass no matter what
      if (junction.data.numpos < minpos_) { continue; }
      // otherwise passed in general vs for denovo
      const bool minreads_passed = junction.data.numreads >= minreads_;
      const bool mindenovo_passed = junction.data.numreads >= mindenovo_;
      if (mindenovo_passed || minreads_passed) {
        // need to update counts_
        size_t contig_idx = contigs_->add(junction.contig.get());
        if (contig_idx >= counts_.size()) {
          // make sure index available for contig
          counts_.resize(contig_idx + 1);
        }
        // get reference to value for contig/coordinates/strand to add
        // make two copies of unstranded junctions to handle experiments with
        // mixed strandedness (unstranded/stranded)
        std::vector<GeneStrandness> strands;
        if (junction.strand == GeneStrandness::AMBIGUOUS) {
          // both strands in sorted order
          if constexpr(GeneStrandness::FORWARD < GeneStrandness::REVERSE) {
            strands = {GeneStrandness::FORWARD, GeneStrandness::REVERSE};
          } else {
            strands = {GeneStrandness::REVERSE, GeneStrandness::FORWARD};
          }
        } else {
          strands = {junction.strand};
        }
        for (GeneStrandness strand : strands) {
          auto& counts = counts_[contig_idx][
            std::make_pair(junction.coordinates, strand)];
          if (minreads_passed) {
            ++(counts.numpassed);
          }
          if (mindenovo_passed) {
            ++(counts.numdenovo);
          }
        }
      }  // done adding junction that passed
    }  // end loop over junctions
    // update num_experiments_
    ++num_experiments_;
  }
};

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONS_HPP
