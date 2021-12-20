/**
 * PassedJunctions.hpp
 *
 * helper classes over GeneJunctions to assign SJJunctions combining 1+ build
 * groups with 1+ experiments
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_PASSEDJUNCTIONS_HPP
#define MAJIQ_PASSEDJUNCTIONS_HPP

#include <vector>
#include <map>
#include <set>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <numeric>
#include <algorithm>
#include <utility>
#include <stdexcept>

#include "MajiqTypes.hpp"
#include "MajiqConstants.hpp"
#include "GeneJunctions.hpp"
#include "SJBinsReads.hpp"
#include "Exons.hpp"
#include "MinExperiments.hpp"


namespace majiq {

namespace detail {
inline bool CloseToPrecedingAnnotatedExon(
    const Exons& exons, const KnownGene& gene, position_t x) {
  // get first exon past x
  auto it = exons.overlap_upper_bound(gene, x);
  if (it == exons.begin()) { return false; }  // no exon overlapping or behind
  do {
    --it;  // decrement to exons behind x (or maybe overlapping first time)
    if (it->gene != gene) {
      // got to next gene, so no exon
      return false;
    } else if (it->coordinates.last_pos() + MAX_DENOVO_DIFFERENCE < x) {
      // we are too far away now
      return false;
    } else if (it->is_full_exon()
        && !(it->is_denovo())
        && it->annotated_coordinates().start <= x
        && x <= it->annotated_coordinates().end + MAX_DENOVO_DIFFERENCE) {
      // we have a close enough annotated exon!
      return true;
    }
  } while (it != exons.begin());
  return false;
}
inline bool CloseToFollowingAnnotatedExon(
    const Exons& exons, const KnownGene& gene, position_t x) {
  // get first exon overlapping or past x
  auto it = exons.overlap_lower_bound(gene, x);
  // keep going forwards until full exon, too far away, or different gene
  for (; it != exons.end(); ++it) {
    if (it->gene != gene) {
      // got to next gene, so no exon
      return false;
    } else if (it->coordinates.first_pos() - MAX_DENOVO_DIFFERENCE > x) {
      // we are too far away now
      return false;
    } else if (it->is_full_exon()
        && !(it->is_denovo())
        && it->annotated_coordinates().start - MAX_DENOVO_DIFFERENCE <= x
        && x <= it->annotated_coordinates().end) {
      // we have a close enough annotated exon!
      return true;
    }
  }
  return false;  // no more exons for any gene
}
}  // namespace detail

class GroupJunctionsGenerator {
 private:
  size_t num_experiments_;
  // known junctions
  const std::shared_ptr<GeneJunctions> known_;
  const std::shared_ptr<Exons> exons_;
  std::vector<size_t> known_num_passed_;
  // unknown denovo junctions per overgene
  std::vector<std::map<GeneJunction, size_t>> contig_denovos_num_passed_;
  // mutexes for exclusive/shared access
  std::shared_mutex group_mutex_;  // shared vs unique access
  std::mutex num_experiments_mutex_;  // num experiments
  // NOTE: we generally expect group_mutex_ to be held befor holding contig_mutex_
  std::unique_ptr<std::mutex[]> contig_mutex_;  // contigs for known/denovo

 public:
  GroupJunctionsGenerator(
      const std::shared_ptr<GeneJunctions>& known,
      const std::shared_ptr<Exons>& exons)
      : num_experiments_{0},
        known_{known},
        exons_{exons},
        known_num_passed_(known_ == nullptr ? 0 : known_->size(), 0),
        contig_denovos_num_passed_(
            known_ == nullptr ? 0 : known_->num_contigs()),
        group_mutex_{},
        num_experiments_mutex_{},
        contig_mutex_{known_ == nullptr ? nullptr
            : std::make_unique<std::mutex[]>(known_->num_contigs())} {
    if (known_ == nullptr) {
      throw std::invalid_argument(
          "GroupJunctionsGenerator requires non-null known GeneJunctions");
    } else if (exons == nullptr) {
      throw std::invalid_argument(
          "GroupJunctionsGenerator requires non-null Exons");
    } else if (known_->parents() != exons_->parents()) {
      throw std::invalid_argument(
          "GroupJunctionsGenerator exons/junctions do not share genes");
    }
  }
  size_t num_experiments() const noexcept { return num_experiments_; }
  size_t num_annotated() const noexcept { return known_num_passed_.size(); }
  size_t num_denovo() const noexcept {
    return std::accumulate(
        contig_denovos_num_passed_.begin(), contig_denovos_num_passed_.end(),
        size_t{0},
        [](size_t s, const auto& x) { return s + x.size(); });
  }
  size_t size() const noexcept { return num_annotated() + num_denovo(); }

  const std::shared_ptr<GeneJunctions>& known() const { return known_; }

  void AddExperiment(const SJJunctionsBins& sjp,
      const ExperimentThresholds& thresholds, bool process_denovo);
  /**
   * Update known junctions in place for passing
   */
  std::shared_ptr<GeneJunctions> UpdateKnownInplace(real_t min_experiments_f) {
    // NOTE: holding group_mutex_ should be equivalent to holding each element
    // of contig_mutex_
    std::scoped_lock lock(group_mutex_, num_experiments_mutex_);
    // determine number of experiments required to pass
    size_t min_experiments = detail::min_experiments_from_float(
        num_experiments_, min_experiments_f);
    for (size_t idx = 0; idx < known_num_passed_.size(); ++idx) {
      if (known_num_passed_[idx] >= min_experiments) {
        (*known_)[idx].passed_build() = true;
      }
    }
    return known_;
  }
  friend class PassedJunctionsGenerator;
};


class PassedJunctionsGenerator {
 private:
  const std::shared_ptr<GeneJunctions> known_;
  std::vector<bool> known_passed_build_;
  std::set<GeneJunction> denovos_passed_build_;
  std::mutex passed_mutex_;

 public:
  size_t num_annotated() const noexcept { return known_passed_build_.size(); }
  size_t num_denovo() const noexcept { return denovos_passed_build_.size(); }
  size_t size() const noexcept { return num_annotated() + num_denovo(); }
  explicit PassedJunctionsGenerator(const std::shared_ptr<GeneJunctions>& known)
      : known_{known},
        known_passed_build_(known_ == nullptr ? 0 : known_->size(), false),
        denovos_passed_build_{},
        passed_mutex_{} {
     if (known_ == nullptr) {
      throw std::invalid_argument(
          "PassedJunctionsGenerator requires non-null known GeneJunctions");
     }
  }

  const std::shared_ptr<GeneJunctions>& known() const { return known_; }

  /**
   * Add junction with gene/coordinates as passed
   */
  void AddJunction(const KnownGene& gene, const OpenInterval& coordinates) {
    // acquire appropriate locks on self
    std::scoped_lock lock(passed_mutex_);
    // junction that we want to add
    GeneJunction j{gene, coordinates, true, true, false};  // passed and denovo
    // check if it is among known junctions
    auto known_it = known_->find(j);
    if (known_it != known_->end()) {
      // it is a known junction, make sure it is marked as passed
      known_passed_build_[known_it - known_->begin()] = true;
    } else {
      // make sure it is among denovos passing build
      denovos_passed_build_.insert(j);
    }
    return;
  }
  void AddJunction(size_t gene_idx, position_t start, position_t end) {
    return AddJunction(
        KnownGene{gene_idx, known_->parents()}, OpenInterval{start, end});
  }

  void AddGroup(
      GroupJunctionsGenerator& group, real_t min_experiments_f) {
    if (known_ != group.known_) {
      throw std::invalid_argument(
          "Added group junctions have different set of known gene junctions");
    }

    // acquire appropriate locks on self and on group junctions generator
    std::scoped_lock lock(
        group.group_mutex_, group.num_experiments_mutex_, passed_mutex_);

    size_t min_experiments = detail::min_experiments_from_float(
        group.num_experiments(), min_experiments_f);

    // update known junctions first
    for (size_t kidx = 0; kidx < known_passed_build_.size(); ++kidx) {
      if (group.known_num_passed_[kidx] >= min_experiments) {
        known_passed_build_[kidx] = true;
      }
    }

    // update denovo junctions
    auto denovo_it = denovos_passed_build_.begin();
    for (const auto& denovos_num_passed : group.contig_denovos_num_passed_) {
      for (const auto& [j, ct] : denovos_num_passed) {
        if (ct < min_experiments) { continue; }  // didn't pass filters
        // get hint for adding in junctions
        denovo_it = std::find_if(denovo_it, denovos_passed_build_.end(),
            [&j](const auto& x) { return j < x; });
        denovos_passed_build_.insert(denovo_it, j);
      }
    }
    return;
  }

  GeneJunctions PassedJunctions(bool denovo_simplified) {
    std::lock_guard lock(passed_mutex_);
    // initialize underlying container of junctions for result
    std::vector<GeneJunction> result_vec(
        known_passed_build_.size() + denovos_passed_build_.size());
    auto r_it = result_vec.begin();
    size_t kidx = 0;  // known
    auto d_it = denovos_passed_build_.begin();  // denovos
    while (r_it != result_vec.end()) {
      // add known junctions before denovos
      for (; kidx < known_passed_build_.size()
          && (d_it == denovos_passed_build_.end() || (*known_)[kidx] < *d_it);
          ++kidx, ++r_it) {
        const auto& j = (*known_)[kidx];
        *r_it = GeneJunction{j.gene, j.coordinates, j.denovo(),
          j.passed_build() || known_passed_build_[kidx], j.simplified()};
      }
      // add denovos before known
      for (; d_it != denovos_passed_build_.end()
          && (kidx == known_passed_build_.size() || *d_it < (*known_)[kidx]);
          ++d_it, ++r_it) {
        // set simplified for junction at d_it
        GeneJunction j = *d_it;
        j.simplified() = denovo_simplified;
        *r_it = std::move(j);
      }
    }  // done filling result_vec
    return GeneJunctions{known_->parents(), std::move(result_vec)};
  }
};

}  // namespace majiq

#endif  // MAJIQ_PASSEDJUNCTIONS_HPP
