/**
 * PassedJunctions.hpp
 *
 * Contig regions for aggregating junctions that pass from multiple experiments
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_PASSEDJUNCTIONS_HPP
#define MAJIQ_PASSEDJUNCTIONS_HPP

#include <vector>
#include <map>
#include <memory>
#include <utility>
#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "MajiqTypes.hpp"
#include "ContigRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "SJJunctions.hpp"


namespace majiq {

// group passed (build, denovo)
struct BuildPassed {
  JunctionPassedStatus status_;
  explicit BuildPassed(JunctionPassedStatus status) : status_{status} { }
  BuildPassed() : BuildPassed{JunctionPassedStatus::NOT_PASSED} { }
  BuildPassed(const BuildPassed&) = default;
  BuildPassed(BuildPassed&&) = default;
  BuildPassed& operator=(const BuildPassed&) = default;
  BuildPassed& operator=(BuildPassed&&) = default;
  BuildPassed& operator+=(const JunctionPassedStatus& rhs) {
    status_ = std::max(status_, rhs);
    return *this;
  }
  BuildPassed& operator+=(const BuildPassed& rhs) {
    return operator+=(rhs.status_);
  }
};

struct GroupCounts {
  uint32_t passed_denovo;
  uint32_t passed_build_only;
  GroupCounts() : passed_denovo{}, passed_build_only{} { }
  explicit GroupCounts(JunctionPassedStatus status) {
    switch (status) {
      case JunctionPassedStatus::DENOVO_PASSED:
        passed_denovo = 1;
        passed_build_only = 0;
        break;
      case JunctionPassedStatus::ANNOTATED_PASSED:
        passed_denovo = 0;
        passed_build_only = 1;
        break;
      case JunctionPassedStatus::NOT_PASSED:
        passed_denovo = 0;
        passed_build_only = 0;
        break;
    }
  }
  GroupCounts& operator+=(const JunctionPassedStatus& rhs) {
    switch (rhs) {
      case JunctionPassedStatus::DENOVO_PASSED:
        ++passed_denovo;
        break;
      case JunctionPassedStatus::ANNOTATED_PASSED:
        ++passed_build_only;
        break;
      case JunctionPassedStatus::NOT_PASSED:
        break;
    }
    return *this;
  }
  GroupCounts& operator+=(const GroupCounts& rhs) {
    passed_denovo += rhs.passed_denovo;
    passed_build_only += rhs.passed_build_only;
    return *this;
  }
  JunctionPassedStatus passed(size_t min_experiments) const {
    if (passed_denovo >= min_experiments) {
      return JunctionPassedStatus::DENOVO_PASSED;
    } else if (passed_denovo + passed_build_only >= min_experiments) {
      return JunctionPassedStatus::ANNOTATED_PASSED;
    } else {
      return JunctionPassedStatus::NOT_PASSED;
    }
  }
};

class PassedJunction : public detail::ContigRegion<OpenInterval, BuildPassed> {
  using BaseT = detail::ContigRegion<OpenInterval, BuildPassed>;
 public:
  const JunctionPassedStatus& status() const noexcept { return data.status_; }
  JunctionPassedStatus& status() noexcept { return data.status_; }
  PassedJunction(KnownContig contig, OpenInterval coordinates,
      GeneStrandness strand, BuildPassed status)
      : BaseT{contig, coordinates, strand, status} { }
  PassedJunction(KnownContig contig, OpenInterval coordinates,
      GeneStrandness strand)
      : PassedJunction{contig, coordinates, strand, BuildPassed{}} { }
  PassedJunction()
      : PassedJunction{KnownContig{}, OpenInterval{}, GeneStrandness{}} { }
  PassedJunction(const PassedJunction&) = default;
  PassedJunction(PassedJunction&&) = default;
  PassedJunction& operator=(const PassedJunction&) = default;
  PassedJunction& operator=(PassedJunction&&) = default;
};
using PassedJunctions = detail::Regions<PassedJunction, true>;

class GroupJunctionsGenerator {
 private:
  // types used
  using interval_stranded_t = std::pair<OpenInterval, GeneStrandness>;
  using interval_stranded_cts_t = std::map<interval_stranded_t, GroupCounts>;

  // parameters for an experiment being marked as passed
  const junction_ct_t minreads_;
  const junction_ct_t mindenovo_;  // REQUIRED >= minreads
  const junction_pos_t minpos_;

  // state of junctions seen
  const std::shared_ptr<Contigs> contigs_;
  std::vector<interval_stranded_cts_t> counts_;
  size_t num_experiments_;

  /**
   * Add junction to dst, using dst_it as starting point for hint.
   * Return iterator following where it was inserted/updated.
   */
  typename interval_stranded_cts_t::iterator AddJunction(
      const SJJunction& junction, interval_stranded_cts_t& dst,
      typename interval_stranded_cts_t::iterator dst_it) {
    JunctionPassedStatus status = junction.passed(
        minreads_, mindenovo_, minpos_);
    if (status == JunctionPassedStatus::NOT_PASSED) { return dst_it; }
    // order of strands to add into contig_cts
    std::vector<GeneStrandness> strand_order;
    if (junction.strand == GeneStrandness::AMBIGUOUS) {
      if constexpr(GeneStrandness::FORWARD < GeneStrandness::REVERSE) {
        strand_order = {GeneStrandness::FORWARD, GeneStrandness::REVERSE};
      } else {
        strand_order = {GeneStrandness::REVERSE, GeneStrandness::FORWARD};
      }
    } else {
      strand_order = {junction.strand};
    }
    // add junction to relevant strands, using hint for fast insertion
    for (const auto& strand : strand_order) {
      const auto insert_key = std::make_pair(junction.coordinates, strand);
      auto insert_value = GroupCounts{status};
      // get hint for insertion for amortized constant insert
      dst_it = std::find_if(dst_it, dst.end(),
          [&insert_key](const auto& x) { return insert_key <= x.first; });
      if (dst_it->first == insert_key) {
        dst_it->second += insert_value;
        // get next iterator since next input will be after current dst_it
        ++dst_it;
      } else {
        dst.insert(dst_it, std::make_pair(insert_key, insert_value));
      }
    }
    return dst_it;
  }

 public:
  GroupJunctionsGenerator(const std::shared_ptr<Contigs>& contigs,
      junction_ct_t minreads, junction_ct_t mindenovo, junction_pos_t minpos)
      : minreads_{minreads},
        mindenovo_{std::max(mindenovo, minreads_)},
        minpos_{minpos},
        contigs_{contigs},
        counts_(contigs_ == nullptr ? 0 : contigs_->size()),
        num_experiments_{} {
    if (contigs_ == nullptr) {
      throw std::invalid_argument(
          "GroupJunctionsGenerator requires non-null contigs");
    }
  }

  void AddExperiment(const SJJunctions& sj) {
    // update num_experiments_
    ++num_experiments_;
    // update counts_
    for (size_t contig_idx = 0; contig_idx < counts_.size(); ++contig_idx) {
      // get contig idx from x matching this one
      const auto opt_sj_contig_idx
        = sj.parents_->safe_idx(contigs_->get(contig_idx));
      if (!opt_sj_contig_idx.has_value()) { continue; }
      const size_t& sj_contig_idx = *opt_sj_contig_idx;
      // destination for junctions from this contig
      auto& dst_cts = counts_[contig_idx];
      auto dst_cts_it = dst_cts.begin();
      // process each junction from the contig, adding to its destination...
      for (auto j_it = sj.begin_parent(sj_contig_idx);
          j_it != sj.end_parent(sj_contig_idx); ++j_it) {
        dst_cts_it = AddJunction(*j_it, dst_cts, dst_cts_it);
      }
    }
  }
  size_t size() const {
    return std::accumulate(counts_.begin(), counts_.end(), size_t{},
        [](size_t s, const interval_stranded_cts_t& x) -> size_t {
        return s + x.size();
        });
  }

  friend class PassedJunctionsGenerator;
};

class PassedJunctionsGenerator {
 private:
  using interval_stranded_t = GroupJunctionsGenerator::interval_stranded_t;
  using interval_stranded_pass_t = std::map<interval_stranded_t, BuildPassed>;

  // junctions that have passed
  const std::shared_ptr<Contigs> contigs_;
  std::vector<interval_stranded_pass_t> passed_;

  static size_t min_experiments_from_float(
      size_t num_experiments, float min_experiments_f) {
    // get min_experiments as size_t from input
    if (min_experiments_f < 0) {
      throw std::invalid_argument("min_experiments must be non-negative");
    } else if (min_experiments_f < 1) {
      // less than 1 ~ percentage of number of experiments
      min_experiments_f *= num_experiments;
    }
    // go to next number of experiments, max being group.num_experiments_
    return std::min(
        static_cast<size_t>(std::ceil(min_experiments_f)),
        num_experiments);
  }

 public:
  explicit PassedJunctionsGenerator(const std::shared_ptr<Contigs>& contigs)
      : contigs_{contigs},
        passed_(contigs_ == nullptr ? 0 : contigs_->size()) {
    if (contigs_ == nullptr) {
      throw std::invalid_argument(
          "PassedJunctionsGenerator requires non-null contigs");
    }
  }
  /**
   * Helper to create group junctions generator with shared set of contigs
   */
  GroupJunctionsGenerator StartGroup(junction_ct_t minreads,
      junction_ct_t mindenovo, junction_pos_t minpos) const {
    return GroupJunctionsGenerator{contigs_, minreads, mindenovo, minpos};
  }
  size_t size() const {
    return std::accumulate(passed_.begin(), passed_.end(), size_t{},
        [](size_t s, const interval_stranded_pass_t& x) -> size_t {
        return s + x.size();
        });
  }
  // TODO(jaicher) add operator+= for PassedJunctions, PassedJunctionsGenerator
  // TODO(jaicher) implement FromPassedJunctions using operator+
  PassedJunctions ToPassedJunctions() const {
    // initialize passed junctions
    std::vector<PassedJunction> result_vec(size());
    size_t jidx = 0;
    for (size_t contig_idx = 0; contig_idx < passed_.size(); ++contig_idx) {
      const interval_stranded_pass_t& contig_passed = passed_[contig_idx];
      const KnownContig contig = (*contigs_)[contig_idx];
      for (const auto& [interval_strand, status] : contig_passed) {
        const auto& [coordinates, strand] = interval_strand;
        result_vec[jidx++]
          = PassedJunction{contig, coordinates, strand, status};
      }
    }
    return PassedJunctions{std::move(result_vec)};
  }
  void AddGroup(const GroupJunctionsGenerator& group, float min_experiments_f) {
    // don't process empty data
    if (group.num_experiments_ == 0 || group.counts_.empty()) { return; }
    // get min_experiments for junction to pass
    size_t min_experiments = min_experiments_from_float(
        group.num_experiments_, min_experiments_f);
    // update passed_ one contig at a time
    for (size_t contig_idx = 0; contig_idx < passed_.size(); ++contig_idx) {
      // get contig_idx from group
      const auto opt_group_contig_idx
        = group.contigs_->safe_idx(contigs_->get(contig_idx));
      if (!opt_group_contig_idx.has_value()) { continue; }
      const size_t& group_contig_idx = *opt_group_contig_idx;
      // get map inserting to and iterator over what's in it so far
      auto& contig_passed = passed_[contig_idx];
      auto passed_it = contig_passed.begin();
      // iterate over group junctions
      const auto& contig_group = group.counts_[group_contig_idx];
      for (auto group_it = contig_group.begin(); group_it != contig_group.end();
          ++group_it) {
        const auto& [key, key_counts] = *group_it;
        auto key_passed = BuildPassed{key_counts.passed(min_experiments)};
        if (key_passed.status_ == JunctionPassedStatus::NOT_PASSED) {
          continue;
        }
        // get first iterator past key we'd like to put in/update (for hint)
        passed_it = std::find_if(passed_it, contig_passed.end(),
            [&key](const auto& x) { return key <= x.first; });
        if (passed_it->first == key) {
          passed_it->second += key_passed;
          // get next iterator since next input will be after current
          ++passed_it;
        } else {
          contig_passed.insert(passed_it, std::make_pair(key, key_passed));
        }
      }
    }
  }
};


}  // namespace majiq


#endif  // MAJIQ_PASSEDJUNCTIONS_HPP
