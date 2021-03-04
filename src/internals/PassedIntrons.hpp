/**
 * PassedIntrons.hpp
 *
 * helper classes over GeneIntrons to assign SJIntronsBins to pass introns
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_PASSEDINTRONS_HPP
#define MAJIQ_PASSEDINTRONS_HPP

#include <algorithm>
#include <memory>
#include <vector>
#include <stdexcept>

#include "GeneIntrons.hpp"
#include "SJJunctionsPositions.hpp"
#include "MinExperiments.hpp"
#include "ExperimentThresholds.hpp"

namespace majiq {

class GroupIntronsGenerator {
 private:
  size_t num_experiments_;
  const std::shared_ptr<GeneIntrons> gene_introns_;
  std::vector<size_t> num_passed_;  // num experiments each intron passed for

 public:
  GroupIntronsGenerator(
      const std::shared_ptr<GeneIntrons>& gene_introns)
      : num_experiments_{},
        gene_introns_{gene_introns},
        num_passed_(gene_introns_ == nullptr ? 0 : gene_introns_->size(), 0) {
    if (gene_introns_ == nullptr) {
      throw std::invalid_argument(
          "GroupIntronsGenerator requires non-null GeneIntrons");
    }
  }

  const std::shared_ptr<GeneIntrons>& introns() { return gene_introns_; }
  const std::vector<size_t>& num_passed() { return num_passed_; }
  size_t num_experiments() const noexcept { return num_experiments_; }
  const std::shared_ptr<GeneIntrons>& gene_introns() const {
    return gene_introns_;
  }
  size_t size() const noexcept { return num_passed_.size(); }

  void AddExperiment(
      const SJIntronsBins& sj, const ExperimentThresholds& exp_thresholds) {
    ++num_experiments_;  // update number of experiments
    const std::shared_ptr<Genes>& genes = gene_introns_->parents();
    // generator of intron thresholds
    IntronThresholdsGenerator thresholds_gen
      = exp_thresholds.intron_thresholds_generator(sj.total_bins());
    // boolean flags of if each gene intron passed, update by iterating over sj
    std::vector<bool> gi_passed(num_passed_.size(), false);
    const ContigIntrons& contig_introns = *(sj.regions());
    for (size_t dst_contig_idx = 0;
        dst_contig_idx < genes->contigs()->size(); ++dst_contig_idx) {
      // try to get matching contigs in sj (since may not share)
      const auto opt_sj_contig_idx = sj.regions()->parents()->safe_idx(
          genes->contigs()->get(dst_contig_idx));
      if (!opt_sj_contig_idx.has_value()) { continue; }
      const size_t& sj_contig_idx = *opt_sj_contig_idx;
      // get iterator to first gene for contig
      KnownGene g_it_start = genes->begin_contig(dst_contig_idx);
      const KnownGene g_it_end = genes->end_contig(dst_contig_idx);
      // for each sj intron on this contig
      for (auto sj_it = contig_introns.begin_parent(sj_contig_idx);
          sj_it != contig_introns.end_parent(sj_contig_idx);
          ++sj_it) {
        // determine if passed
        constexpr junction_pos_t NO_STACKS = 0;  // TODO(jaicher) check stacks?
        if (!sj.passed(
              sj_it - contig_introns.begin(), thresholds_gen, NO_STACKS)) {
          continue;
        }
        // get first gene that doesn't precede current intron
        g_it_start = std::find_if(g_it_start, g_it_end,
            [&sj_coords = sj_it->coordinates](const KnownGene& x) {
            return !IntervalPrecedes(x.coordinates(), sj_coords); });
        if (g_it_start == g_it_end) { break; }  // no genes on contig remain
        // then loop over genes starting from there until intron precedes genes
        for (KnownGene g_it = g_it_start;
            g_it != g_it_end
            && !IntervalPrecedes(sj_it->coordinates, g_it.coordinates());
            ++g_it) {
          if (
              // strand of contig intron matches gene
              (sj_it->strand == GeneStrandness::AMBIGUOUS
               || sj_it->strand == g_it.strand())
              // and intron is contained by gene
              && IntervalSubsets(sj_it->coordinates, g_it.coordinates())) {
            // find overlapping gene introns and mark them as passed
            for (
                auto it = gene_introns_->overlap_lower_bound(
                  // overlap starting 1 before to handle zero-length introns
                  // but means that first might not actually intersect
                  g_it, sj_it->coordinates.start - 1);
                it != gene_introns_->end_parent(g_it)
                && !IntervalPrecedes(sj_it->coordinates, it->coordinates);
                ++it) {
              if (
                  IntervalIntersects(sj_it->coordinates, it->coordinates)
                  // if contig intron is annotated, so must the gene intron
                  && (!sj_it->annotated() || !it->denovo())) {
                // it overlaps a passed intron --> mark as passed
                gi_passed[it - gene_introns_->begin()] = true;
              }
            }  // loop over intersecting introns in a gene
          }  // only if strand matches and intron contained by gene
        }  // for all genes that don't come after intron
      }  // for each sj intron on the contig
    }  // for each contig
    // update num_passed_
    for (size_t i = 0; i < num_passed_.size(); ++i) {
      if (gi_passed[i]) { ++num_passed_[i]; }
    }
  }

  // updates GeneIntrons in place, resets num_passed and num_experiments
  const std::shared_ptr<GeneIntrons>& UpdateInplace(real_t min_experiments_f) {
    size_t min_experiments = detail::min_experiments_from_float(
        num_experiments_, min_experiments_f);
    // update introns that had enough experiments to pass
    for (size_t idx = 0; idx < num_passed_.size(); ++idx) {
      if (num_passed_[idx] >= min_experiments) {
        (*gene_introns_)[idx].passed_build() = true;
      }
    }
    // clear num_experiments, num_passed
    num_experiments_ = 0;
    std::fill(num_passed_.begin(), num_passed_.end(), 0);
    // return pointer to updated introns
    return gene_introns_;
  }
};

}  // namespace majiq

#endif  // MAJIQ_PASSEDINTRONS_HPP
