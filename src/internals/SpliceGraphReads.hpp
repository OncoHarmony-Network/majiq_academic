/**
 * SpliceGraphReads.hpp
 *
 * Get number of reads for splicegraph junctions and introns from an
 * experiments' SJ junctions/introns
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_SPLICEGRAPHREADS_HPP
#define MAJIQ_SPLICEGRAPHREADS_HPP

#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "GeneIntrons.hpp"
#include "GeneJunctions.hpp"
#include "SJBinsReads.hpp"


namespace majiq {

class SpliceGraphReads {
 private:
  const std::shared_ptr<GeneIntrons> introns_;
  const std::shared_ptr<GeneJunctions> junctions_;
  const std::vector<real_t> introns_reads_;
  const std::vector<real_t> junctions_reads_;

 public:
  SpliceGraphReads(
      const std::shared_ptr<GeneIntrons>& introns,
      const std::shared_ptr<GeneJunctions>& junctions,
      std::vector<real_t>&& introns_reads,
      std::vector<real_t>&& junctions_reads)
      : introns_{introns},
        junctions_{junctions},
        introns_reads_{std::move(introns_reads)},
        junctions_reads_{std::move(junctions_reads)} {
    if (introns_ == nullptr) {
      throw std::runtime_error("SpliceGraphReads given null introns");
    } else if (junctions_ == nullptr) {
      throw std::runtime_error("SpliceGraphReads given null junctions");
    } else if (introns_->size() != introns_reads_.size()) {
      throw std::runtime_error(
          "SpliceGraphReads introns reads do not match introns in size");
    } else if (junctions_->size() != junctions_reads_.size()) {
      throw std::runtime_error(
          "SpliceGraphReads junctions reads do not match junctions in size");
    }
  }
  SpliceGraphReads(const SpliceGraphReads&) = default;
  SpliceGraphReads(SpliceGraphReads&&) = default;
  SpliceGraphReads& operator=(const SpliceGraphReads&) = delete;
  SpliceGraphReads& operator=(SpliceGraphReads&&) = delete;

  const std::shared_ptr<GeneIntrons>& introns() const { return introns_; }
  const std::shared_ptr<GeneJunctions>& junctions() const { return junctions_; }
  const std::vector<real_t>& introns_reads() const { return introns_reads_; }
  const std::vector<real_t>& junctions_reads() const {
    return junctions_reads_;
  }

  // how do we get this from SJ files?
 private:
  template <typename GeneRegionsT, typename SJBinsT>
  static std::vector<real_t> annotate_numreads(
      const GeneRegionsT& regions, const SJBinsT& sj) {
    static_assert(
        (std::is_same_v<GeneRegionsT, GeneIntrons>
         && std::is_same_v<SJBinsT, SJIntronsBins>)
        || (std::is_same_v<GeneRegionsT, GeneJunctions>
          && std::is_same_v<SJBinsT, SJJunctionsBins>),
        "annotate_numreads needs matching gene regions and SJ bins");
    constexpr bool IS_INTRON = std::is_same_v<SJIntronsBins, SJBinsT>;
    using RegionIntervalBeforeT = std::conditional_t<
      IS_INTRON,
      // zero-length introns [a, a-1], [b+1, b] overlap [a, b]
      IntervalPrecedesT<true>,
      std::less<std::conditional_t<IS_INTRON, ClosedInterval, OpenInterval>>>;

    const auto& sj_regions = *(sj.regions());
    std::vector<junction_pos_t> weights(regions.size(), 0);
    std::vector<real_t> result(regions.size(), 0);

    // contigs (which might not be the same between gene regions/sj regions
    Genes& genes = *(regions.parents());
    const Contigs& genes_contigs = *(genes.parents());
    const Contigs& sj_contigs = *(sj_regions.parents());
    // loop over contigs with same id
    for (size_t dst_contig_idx = 0;
        dst_contig_idx < genes_contigs.size(); ++dst_contig_idx) {
      const auto opt_sj_contig_idx = sj_contigs.safe_idx(
          genes_contigs.get(dst_contig_idx));
      if (!opt_sj_contig_idx.has_value()) { continue; }
      const size_t& sj_contig_idx = *opt_sj_contig_idx;

      // get iterator over sj regions for contig
      auto sj_it = sj_regions.begin_parent(sj_contig_idx);
      const auto sj_it_end = sj_regions.end_parent(sj_contig_idx);
      // no sj regions for contig?
      if (sj_it == sj_it_end) { continue; }

      // get indexes over regions in contig sorted order for this contig
      std::vector<size_t> region_idx;
      {
        auto region_it = regions.begin_parent(genes.begin_parent(dst_contig_idx));
        const auto region_it_end
          = regions.begin_parent(genes.end_parent(dst_contig_idx));
        // no gene regions for contig?
        if (region_it == region_it_end) { continue; }
        // set up indexes into regions for the contig, then sort by coordinates
        region_idx.resize(region_it_end - region_it);
        std::iota(
            region_idx.begin(), region_idx.end(), region_it - regions.begin());
        std::sort(
            region_idx.begin(), region_idx.end(),
            [&x = regions](size_t i, size_t j) {
            return x[i].coordinates < x[j].coordinates; });
      }

      // loop over indexes vs sj to find overlaps and update weights, counts
      auto idx_it = region_idx.begin();
      const auto idx_it_end = region_idx.end();
      for (; sj_it != sj_it_end; ++sj_it) {
        // update idx_it to first index that could overlap given sorting
        idx_it = std::find_if(idx_it, idx_it_end,
            [&x = regions, &sj_coordinates = sj_it->coordinates](size_t i) {
            return !RegionIntervalBeforeT{}(
                x[i].coordinates, sj_coordinates); });
        if (idx_it == idx_it_end) { break; }  // no more regions, do next contig

        // only calculate weight/numreads once, on first overlap
        junction_pos_t weight = 0;
        real_t numreads = 0;
        // find overlaps
        for (auto overlap_it = idx_it;
            overlap_it != idx_it_end
            && !RegionIntervalBeforeT{}(
              sj_it->coordinates, regions[*overlap_it].coordinates);
            ++overlap_it) {
          // potential match
          const auto& region = regions[*overlap_it];
          // ignore non-matching strand
          if (sj_it->strand != GeneStrandness::AMBIGUOUS
              && sj_it->strand != region.strand()) { continue; }
          // calculate weight/numreads if first time
          if (weight == 0) {
            const auto sj_idx = sj_it - sj_regions.begin();
            weight = sj.weight(sj_idx);
            constexpr junction_pos_t NUM_STACKS = 0;  // raw counts: no stacks
            numreads = sj.numreads(sj_idx, NUM_STACKS);
          }
          weights[*overlap_it] += weight;
          result[*overlap_it] += numreads;
        }  // loop over overlapping regions for given sj
      }  // loop over sj in the contig
    }  // loop over contigs

    // normalize numreads
    const real_t total_bins_f = static_cast<real_t>(sj.total_bins());
    for (size_t i = 0; i < weights.size(); ++i) {
      const auto& weight = weights[i];
      if (weight == 0) { continue; }
      result[i] *= total_bins_f / weight;
    }  // normalize numreads
    return result;
  }

 public:
  static SpliceGraphReads FromSJ(
      const std::shared_ptr<GeneIntrons>& introns,
      const std::shared_ptr<GeneJunctions>& junctions,
      const SJIntronsBins& sj_introns,
      const SJJunctionsBins& sj_junctions) {
    std::vector<real_t> introns_reads
      = annotate_numreads(*introns, sj_introns);
    std::vector<real_t> junctions_reads
      = annotate_numreads(*junctions, sj_junctions);
    return SpliceGraphReads{
      introns, junctions, std::move(introns_reads), std::move(junctions_reads)};
  }
};

}  // namespace majiq


#endif  // MAJIQ_SPLICEGRAPHREADS_HPP
