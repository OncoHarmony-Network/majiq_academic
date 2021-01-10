/**
 * OverGenes.hpp
 *
 * Identify offsets and intervals for contigs, overlapping genes
 *
 * Copyright 2020 <University of Pennsylvania
 */
#ifndef MAJIQ_OVERGENES_HPP
#define MAJIQ_OVERGENES_HPP

#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <utility>

#include "Interval.hpp"
#include "Regions.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
using OverGene = detail::ContigRegion<ClosedInterval>;

class OverGenes {
 public:
  using VecOverGene = std::vector<OverGenes>;

 private:
  std::vector<size_t> contig_offsets_;  // into overgenes (length 1+ contigs)
  std::vector<OverGene> overgenes_;
  std::vector<size_t> overgene_offsets_;  // into genes (length 1+ overgenes)

 public:
  inline size_t num_contigs() const noexcept {
    return contig_offsets_.size() - 1;
  }
  inline size_t contig_overgenes(size_t contig_idx) const noexcept {
    return contig_idx < num_contigs()
      ? contig_offsets_[contig_idx + 1] - contig_offsets_[contig_idx]
      : 0;
  }
  inline std::pair<size_t, size_t> overgenes_range(size_t contig_idx) const {
    return std::make_pair(
        contig_offsets_[contig_idx], contig_offsets_[contig_idx + 1]);
  }
  inline size_t size() const noexcept { return overgenes_.size(); }
  inline const OverGene& get_overgene(size_t overgene_idx) const {
    return overgenes_[overgene_idx];
  }
  inline std::pair<size_t, size_t> genes_range(size_t overgene_idx) const {
    return std::make_pair(
        overgene_offsets_[overgene_idx], overgene_offsets_[overgene_idx + 1]);
  }

  // access for pybind
  const std::vector<size_t>& contig_offsets() { return contig_offsets_; }
  const std::vector<OverGene>& overgenes() { return overgenes_; }
  const std::vector<size_t>& overgene_offsets() { return overgene_offsets_; }

  /**
   * overgene_idx identifies the index of the unique overgene that contains
   * the specified region, or -1 if no such region exists
   */
  template <class IntervalT>
  size_t overgene_idx(
      size_t contig_idx, const IntervalT& x) const noexcept {
    static_assert(std::is_base_of<detail::Interval, IntervalT>::value,
        "overgene_idx requires IntervalT derived from detail::Interval");
    if (contig_idx >= num_contigs()) {
      return -1;
    }
    // iterators for range of overgenes for the contig
    auto contig_begin = overgenes_.begin() + contig_offsets_[contig_idx];
    auto contig_end = overgenes_.begin() + contig_offsets_[contig_idx + 1];
    // iterator to first interval that is greater than x in the range
    auto ub = std::upper_bound(contig_begin,  contig_end, x,
        [](const OverGene& lhs, const detail::Interval& rhs) -> bool {
          return lhs.coordinates < rhs;
        });
    // x < *ub, so *ub overlaps if and only if they share a start coordinate.
    // Otherwise, *(ub -1) is the only candidate
    if (ub < contig_end && IntervalSubsets(x, *ub)) {
      return ub - overgenes_.begin();
    } else if (ub > contig_begin && IntervalSubsets(x, *(ub - 1))) {
      return (ub - 1) - overgenes_.begin();
    } else {
      return -1;
    }
  }
  template <class IntervalT>
  size_t overgene_idx(const detail::ContigRegion<IntervalT>& x) const noexcept {
    return overgene_idx(x.contig.contig_idx, x.coordinates);
  }
  template <class IntervalT>
  size_t overgene_idx(const detail::GeneRegion<IntervalT>& x) const noexcept {
    return overgene_idx(x.gene.get().contig_idx, x.coordinates);
  }

  // constructors
  OverGenes(
      const std::shared_ptr<Contigs>& contigs,
      const std::shared_ptr<Genes>& genes)
      : contig_offsets_{1},
        overgenes_{},
        overgene_offsets_{1} {
    if (!genes->is_sorted()) {
      throw std::runtime_error("OverGenes requires genes in sorted order");
    }
    // current contig
    size_t contig_idx = 0;

    // current overgene
    // overgene_idx = overgene.size()
    constexpr position_t EMPTY = -1;
    position_t og_start = EMPTY;
    position_t og_end = EMPTY;

    for (size_t gene_idx = 0; gene_idx < genes->size(); ++gene_idx) {
      const Gene& gene = genes->get(gene_idx);
      if (
          (contig_idx < gene.contig.contig_idx || og_end < gene.interval.start)
          && og_start != EMPTY
          && og_end != EMPTY) {
        // add an overgene
        overgenes_.push_back(
            OverGene{(*contigs)[contig_idx], ClosedInterval{og_start, og_end},
              GeneStrandness::AMBIGUOUS});
        overgene_offsets_.push_back(gene_idx);
        // reset og_start, og_end
        og_start = gene.interval.start;
        og_end = gene.interval.end;
      } else {
        // start shouldn't change based on ordering, but og_end could
        og_end = std::max(og_end, gene.interval.end);
      }
      while (contig_idx < gene.contig.contig_idx) {
        // add/note end of contig
        contig_offsets_.push_back(overgenes_.size());
        ++contig_idx;
      }
    }
    // add last overgene, if any
    if (og_start != EMPTY && og_end != EMPTY) {
      overgenes_.push_back(
          OverGene{(*contigs)[contig_idx], ClosedInterval{og_start, og_end},
            GeneStrandness::AMBIGUOUS});
      overgene_offsets_.push_back(genes->size());
    }
    // fill in any remaining contigs
    while (contig_idx < contigs->size()) {
      // add/note end of contig
      contig_offsets_.push_back(overgenes_.size());
      ++contig_idx;
    }
    return;
  }
  OverGenes(const OverGenes& x) = default;
  OverGenes(OverGenes&& x) = default;
  OverGenes& operator=(const OverGenes& x) = default;
  OverGenes& operator=(OverGenes&& x) = default;
};
}  // namespace majiq

#endif  // MAJIQ_OVERGENES_HPP
