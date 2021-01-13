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
#include <sstream>
#include <stdexcept>

#include "MajiqTypes.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"

namespace majiq {
// identity/read count for nonzero junction positions
struct PositionReads {
  junction_pos_t pos;
  junction_ct_t reads;
};
inline bool operator<(const PositionReads& x, const PositionReads& y) noexcept {
  return x.pos < y.pos;
}

// total reads and nonzero positions for a junction
struct ExperimentCounts {
  junction_ct_t numreads;
  junction_pos_t numpos;
};

// number of experiments passing filters
struct GroupCounts {
  uint32_t numpassed;  // passing build filters
  uint32_t numdenovo;  // passing denovo filters
};

using ContigJunction = detail::ContigRegion<OpenInterval, detail::EmptyDataT>;
using SJJunction = detail::ContigRegion<OpenInterval, ExperimentCounts>;
using GroupJunction = detail::ContigRegion<OpenInterval, GroupCounts>;

using SJJunctions = detail::ContigRegions<SJJunction>;
using GroupJunctions = detail::ContigRegion<GroupJunction>;

// how do we get from SJJunctions to GroupJunctions?
class GroupJunctionsGenerator {
 private:
  const junction_ct_t minreads_;
  const junction_ct_t mindenovo_;  // REQUIRED >= minreads
  const junction_pos_t minpos_;

  std::shared_ptr<Contigs> contigs_;

  size_t num_experiments_;
  std::vector<  // over contigs
    std::map<std::pair<OpenInterval, GeneStrandness>, GroupCounts>> counts_;

 public:
  GroupJunctionsGenerator(junction_ct_t minreads, junction_ct_t mindenovo,
      junction_pos_t minpos, std::shared_ptr<Contigs> contigs)
      : minreads_{minreads}, mindenovo_{mindenovo}, minpos_{minpos},
        contigs_{contigs == nullptr ? std::make_shared<Contigs>() : contigs} {
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

  // filter junctions that passed group filters
  GroupJunctions PassedJunctions(size_t min_experiments) {
    std::vector<GroupJunction> passed;
    for (size_t contig_idx = 0; contig_idx < contigs_->size(); ++contig_idx) {
      for (auto&& [coord_strand, ct] : counts_[contig_idx]) {
        if (ct.numpassed >= min_experiments) {
          passed.emplace_back(
              (*contigs_)[contig_idx], coord_strand.first, coord_strand.second,
              ct);
        }
      }
    }
    // create GroupJunctions
    return GroupJunctions{passed};
  }
  // add experiment
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
        // get reference to value for contig/coordinates/strand
        auto& counts = counts_[contig_idx][
          std::make_pair(junction.coordinates, junction.strand)];
        if (minreads_passed) {
          ++(counts.numpassed);
        }
        if (mindenovo_passed) {
          ++(counts.numdenovo);
        }
      }  // done adding junction that passed
    }  // end loop over junctions
    // update num_experiments_
    ++num_experiments_;
  }
};

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONS_HPP
