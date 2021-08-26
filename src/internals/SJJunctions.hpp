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
#include <algorithm>

#include "MajiqTypes.hpp"
#include "ContigRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "ExperimentThresholds.hpp"

namespace majiq {

class SJJunction : public detail::ContigRegion<OpenInterval, EmptyDataT> {
  using BaseT = detail::ContigRegion<OpenInterval, EmptyDataT>;

 public:
  SJJunction(KnownContig contig, OpenInterval coordinates,
      GeneStrandness strand)
      : BaseT{contig, coordinates, strand} { }
  SJJunction() : SJJunction{KnownContig{}, OpenInterval{}, GeneStrandness{}} { }
  SJJunction(const SJJunction&) = default;
  SJJunction(SJJunction&&) = default;
  SJJunction& operator=(const SJJunction&) = default;
  SJJunction& operator=(SJJunction&&) = default;
};

class SJJunctions : public detail::Regions<SJJunction, true> {
  using BaseT = detail::Regions<SJJunction, true>;

 public:
  SJJunctions(
      const std::shared_ptr<Contigs>& contigs, std::vector<SJJunction>&& x)
      : BaseT{contigs, std::move(x)} {
    if (parents() == nullptr) {
      throw std::invalid_argument("SJJunctions cannot have null contigs");
    }
  }

  SJJunctions ToUnstranded() const {
    // accumulate unstranded junctions
    std::vector<SJJunction> unstranded;
    // loop through junctions in contig stranded order
    for (auto base_it = begin(); base_it != end();) {
      // get next junction that's different when ignoring strand
      auto next_it = std::find_if(std::next(base_it), end(),
          [&base = *base_it](const SJJunction& x) {
          // first junction where base < x (since ignoring strand)
          return detail::CompareContigUnstranded<SJJunction>()(base, x); });
      // push current base junction into unstranded, set strand to ambiguous
      unstranded.emplace_back(
          base_it->contig, base_it->coordinates, GeneStrandness::AMBIGUOUS);
      // update base_it to next_it
      base_it = next_it;
    }
    // construct new SJJunctions using unstranded junctions we have accumulated
    return SJJunctions(parents(), std::move(unstranded));
  }

  SJJunctions FlipStrand() const {
    // copy junctions from self, flipping strand
    std::vector<SJJunction> sj_vec;
    sj_vec.reserve(size());
    std::transform(begin(), end(), std::back_inserter(sj_vec),
        [](const SJJunction& x) {
        return SJJunction{x.contig, x.coordinates, FlipGeneStrand(x.strand)}; });
    // this may not be in contig stranded order (because strands are flipped)
    // we know that it's in contig unstranded order, so we iterate to find
    // intervals [base_it, next_it) that are equivalent other than strand, then
    // sort by strand
    for (auto base_it = sj_vec.begin(); base_it != sj_vec.end();) {
      // get next junction that's different when ignoring strand
      auto next_it = std::find_if(std::next(base_it), sj_vec.end(),
          [&base = *base_it](const SJJunction& x) {
          // first junction where base < x (since ignoring strand)
          return detail::CompareContigUnstranded<SJJunction>()(base, x); });
      // sort [base_it, next_it) by strand
      std::sort(base_it, next_it, [](const SJJunction& x, const SJJunction& y) {
          return x.strand < y.strand; });
      // update base_it to next_it
      base_it = next_it;
    }
    // construct new SJJunctions using unstranded junctions we have accumulated
    return SJJunctions(parents(), std::move(sj_vec));
  }
};

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONS_HPP
