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
/**
 * Given position x, determine if gene has close exon behind (or overlapping)
 */
inline bool CloseToPrecedingAnnotatedExon(
    const Exons& exons, const KnownGene& gene, position_t x);

/**
 * Given position x, determine if gene has close exon after (or overlapping)
 */
inline bool CloseToFollowingAnnotatedExon(
    const Exons& exons, const KnownGene& gene, position_t x);

/**
 * Assign denovo junction junction to potential matched genes in range [first,
 * last)
 *
 * Track genes in og that junction could be assigned to. So, we prioritize
 * junctions that are close (or overlapping) an exon in the gene, preferring
 * genes where the junction is close on both ends over just one end vs neither
 */
inline void AssignDenovoJunction(const SJJunction& junction, const KnownGene& first,
    const KnownGene& last, const Exons& exons,
    std::map<GeneJunction, size_t>& dst);
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

  /**
   * Add junction bin reads as experiment to GroupJunctionsGenerator
   */
  inline void AddExperiment(const SJJunctionsBins& sjp,
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

// implementation of detail::CloseToPrecedingAnnotatedExon
inline bool detail::CloseToPrecedingAnnotatedExon(
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

// implementation of detail::CloseToFollowingAnnotatedExon
inline bool detail::CloseToFollowingAnnotatedExon(
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

// implementation of detail::AssignDenovoJunction
void detail::AssignDenovoJunction(
    const SJJunction& junction, const KnownGene& first,
    const KnownGene& last, const Exons& exons,
    std::map<GeneJunction, size_t>& dst) {
  // when junction is close to exons on both junction ends as appropriate
  std::vector<KnownGene> close_x2;
  // when only one end of junction close to a gene exon
  std::vector<KnownGene> close_x1;
  // when within gene boundaries but not close to exons at all
  std::vector<KnownGene> far;

  // check all of the overgenes
  for (KnownGene gene = first; gene != last; ++gene) {
    // if strand doesn't match, ignore
    if (!(junction.strand == GeneStrandness::AMBIGUOUS
          || junction.strand == gene.strand())) {
      continue;
    }
    // if junction cannot be placed inside of gene, ignore
    {
      const auto exon_first = exons.begin_parent(gene);
      const auto exon_last = exons.end_parent(gene);
      if (exon_first == exon_last
          || !IntervalSubsets(junction.coordinates,
            ClosedInterval{exon_first->coordinates.first_pos(),
              std::prev(exon_last)->coordinates.last_pos()})) {
        continue;
      }
    }
    // count how many of start/end are close to exon
    int n_close = 0;  // how many of start/end close to an exon?
    if (detail::CloseToPrecedingAnnotatedExon(
          exons, gene, junction.coordinates.start)) {
      ++n_close;
    }
    if (detail::CloseToFollowingAnnotatedExon(
          exons, gene, junction.coordinates.end)) {
      ++n_close;
    }
    // priority of this gene/where we track it
    auto& gene_priority_vec
      = n_close == 0 ? far : (n_close == 1 ? close_x1 : close_x2);
    gene_priority_vec.push_back(gene);
  }  // done looking at all potential genes
  // so get the highest priority vector of genes that is nonempty
  const auto& matched_genes
    = close_x2.empty() ? (close_x1.empty() ? far : close_x1) : close_x2;
  // add count to dst for this junction in those genes, if any
  std::for_each(matched_genes.begin(), matched_genes.end(),
      [&junction, &dst](const KnownGene& gene) {
      auto key = GeneJunction{gene, junction.coordinates, true, true, false};
      ++dst[key];
      });
  return;
}

// implementation of GroupJunctionsGenerator::AddExperiment
void GroupJunctionsGenerator::AddExperiment(const SJJunctionsBins& sjp,
    const ExperimentThresholds& thresholds, bool process_denovo) {
  std::shared_lock group_lock{group_mutex_};
  {  // increment number of experiments
    std::lock_guard num_experiments_lock{num_experiments_mutex_};
    ++num_experiments_;
  }
  const SJJunctions& sj = *(sjp.regions());
  if (sj.empty()) { return; }  // no junctions to add

  const std::shared_ptr<Genes>& genes = exons_->parents();

  // for each contig we have genes for
  for (size_t dst_contig_idx = 0;
      dst_contig_idx < genes->contigs()->size();
      ++dst_contig_idx) {
    // try to get matching contigs in sj (since often does not share contigs)
    const auto opt_sj_contig_idx = sj.parents()->safe_idx(
        genes->contigs()->get(dst_contig_idx));
    if (!opt_sj_contig_idx.has_value()) { continue; }
    const size_t& sj_contig_idx = *opt_sj_contig_idx;

    // get gene junctions (or indexes to them) in contig sorted order
    // so we iterate over sj, gene junctions (known) and genes (for denovo)
    auto sj_it = sj.begin_parent(sj_contig_idx);
    const auto sj_it_end = sj.end_parent(sj_contig_idx);
    if (sj_it == sj_it_end) { continue; }  // no sj junctions for contig
    KnownGene g_it_start = genes->begin_contig(dst_contig_idx);
    const KnownGene g_it_end = genes->end_contig(dst_contig_idx);
    if (g_it_start == g_it_end) { continue; }  // no genes for contig

    // (indexes to) known genejunctions for this contig sorted by coordinates
    std::vector<size_t> known_idx(
        // number of known junctions for the contig
        known_->parent_idx_offsets_[g_it_end.idx_]
        - known_->parent_idx_offsets_[g_it_start.idx_]);
    std::iota(known_idx.begin(), known_idx.end(),
        known_->parent_idx_offsets_[g_it_start.idx_]);
    std::sort(known_idx.begin(), known_idx.end(),
        [&x = *known_](size_t i, size_t j) {
        // ignore gene, strand information, just coordinates within a contig
        return x[i].coordinates < x[j].coordinates; });
    auto k_it_start = known_idx.begin();

    // get denovos num passed for specified contig
    auto& denovos_num_passed = contig_denovos_num_passed_[dst_contig_idx];

    // lock on contig being updated
    std::lock_guard contig_lock{contig_mutex_[dst_contig_idx]};

    for (; sj_it != sj_it_end; ++sj_it) {
      const SJJunction& junction = *sj_it;
      JunctionPassedStatus status = sjp.passed(sj_it - sj.begin(), thresholds);
      if (status == JunctionPassedStatus::NOT_PASSED) { continue; }
      // advance to first known junction that could overlap with sj_it
      k_it_start = std::find_if(k_it_start, known_idx.end(),
          [&x = *known_, &junction](size_t i) {
          return !(x[i].coordinates < junction.coordinates); });
      bool found = false;  // was junction found in known?
      // try to find matching known junction
      for (auto k_it = k_it_start;
          // still have known junctions
          k_it != known_idx.end()
          // that have matching coordinates to junction
          && (*known_)[*k_it].coordinates == junction.coordinates;
          ++k_it) {
        const GeneJunction& kj = (*known_)[*k_it];  // gene junction for k_it
        // if strand matches and junction passed status ~ denovo status
        if ((junction.strand == GeneStrandness::AMBIGUOUS
             || junction.strand == kj.gene.strand())
            && (!kj.denovo() || status == JunctionPassedStatus::DENOVO_PASSED)
            ) {
          ++known_num_passed_[*k_it];
          found = true;
        }
      }  // loop over potentially matching known junctions
      if (!process_denovo) {
        // if no denovos, end of known junctions in contig is stopping condition
        if (k_it_start == known_idx.end()) { break; }
        // otherwise, we are processing junctions not previously known to us
      } else if (!found && status == JunctionPassedStatus::DENOVO_PASSED) {
        // advance to first gene that could overlap with junction
        g_it_start = std::find_if(g_it_start, g_it_end,
            [&junction](const KnownGene& x) {
            return !IntervalPrecedes(x.coordinates(), junction.coordinates); });
        if (g_it_start == g_it_end) { break; }  // no more genes for contig
        // gene can be from there until first gene that starts after junction
        auto g_it_after = std::find_if(g_it_start, g_it_end,
            [&junction](const KnownGene& x) {
            return IntervalPrecedes(junction.coordinates, x.coordinates()); });
        if (g_it_start != g_it_after) {
          detail::AssignDenovoJunction(
              junction, g_it_start, g_it_after, *exons_, denovos_num_passed);
        }
      }
    }  // end loop over sj junctions on contig
  }  // end loop over contig
}

}  // namespace majiq

#endif  // MAJIQ_PASSEDJUNCTIONS_HPP
