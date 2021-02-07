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
#include <stdexcept>

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


using ContigJunction = detail::ContigRegion<OpenInterval, EmptyDataT>;
using SJJunction = detail::ContigRegion<OpenInterval, ExperimentCounts>;
using GroupJunction = detail::ContigRegion<OpenInterval, GroupPassed>;

using SJJunctions = detail::Regions<SJJunction, true>;
using GroupJunctions = detail::Regions<GroupJunction, true>;


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
          const auto& [coordinates, strand] = coord_strand;
          passed.emplace_back(
              (*contigs_)[contig_idx], coordinates, strand,
              GroupPassed{ct.numpassed >= min_experiments,
                          ct.numdenovo >= min_experiments});
        }
      }
    }
    // create GroupJunctions
    return GroupJunctions{std::move(passed)};
  }
  /**
   * Add junctions from experiment to group towards group filters
   */
  void AddExperiment(const SJJunctions& experiment);
};


}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONS_HPP
