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
#include <algorithm>
#include <utility>
#include <stdexcept>

#include "MajiqTypes.hpp"
#include "MajiqConstants.hpp"
#include "GeneJunctions.hpp"
#include "SJJunctions.hpp"
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
  std::map<GeneJunction, size_t> denovos_num_passed_;

 public:
  GroupJunctionsGenerator(
      const std::shared_ptr<GeneJunctions>& known,
      const std::shared_ptr<Exons>& exons)
      : num_experiments_{0},
        known_{known},
        exons_{exons},
        known_num_passed_(known->size(), 0),
        denovos_num_passed_{} {
    if (known_->parents() != exons_->parents()) {
      throw std::invalid_argument(
          "GroupJunctionsGenerator exons/junctions do not share genes");
    }
  }
  size_t num_experiments() const noexcept { return num_experiments_; }
  size_t num_annotated() const noexcept { return known_num_passed_.size(); }
  size_t num_denovo() const noexcept { return denovos_num_passed_.size(); }
  size_t size() const noexcept { return num_annotated() + num_denovo(); }

  void AddExperiment(const SJJunctions& sj, junction_ct_t minreads,
      junction_ct_t mindenovo, junction_pos_t minpos, bool process_denovo);
  /**
   * Update known junctions in place for passing
   */
  std::shared_ptr<GeneJunctions> UpdateKnownInplace(float min_experiments_f) {
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

 public:
  size_t num_annotated() const noexcept { return known_passed_build_.size(); }
  size_t num_denovo() const noexcept { return denovos_passed_build_.size(); }
  size_t size() const noexcept { return num_annotated() + num_denovo(); }
  explicit PassedJunctionsGenerator(const std::shared_ptr<GeneJunctions>& known)
      : known_{known},
        known_passed_build_(known->size(), false),
        denovos_passed_build_{} { }

  void AddGroup(const GroupJunctionsGenerator& group, float min_experiments_f) {
    if (known_ != group.known_) {
      throw std::invalid_argument(
          "Added group junctions have different set of known gene junctions");
    }

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
    for (const auto& [j, ct] : group.denovos_num_passed_) {
      if (ct < min_experiments) { continue; }  // didn't pass filters
      // get hint for adding in junctions
      denovo_it = std::find_if(denovo_it, denovos_passed_build_.end(),
          [&j](const auto& x) { return j < x; });
      denovos_passed_build_.insert(denovo_it, j);
    }
    return;
  }

  GeneJunctions PassedJunctions() const {
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
        *r_it = *d_it;
      }
    }  // done filling result_vec
    return GeneJunctions{known_->parents(), std::move(result_vec)};
  }
};

}  // namespace majiq

#endif  // MAJIQ_PASSEDJUNCTIONS_HPP
