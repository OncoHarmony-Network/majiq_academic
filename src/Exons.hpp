/**
 * Exons.hpp
 *
 * Exons for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_EXONS_HPP
#define MAJIQ_EXONS_HPP

#include <tuple>
#include <functional>
#include <boost/functional/hash.hpp>

#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
struct Exon : detail::GeneRegion<ClosedInterval> {
 public:
  // constructors
  Exon(KnownGene _gene, ClosedInterval _coordinates)
      : detail::GeneRegion<ClosedInterval>{_gene, _coordinates} {
  }
  Exon(const Exon& x) = default;
  Exon(Exon&& x) = default;
  Exon& operator=(const Exon& x) = default;
  Exon& operator=(Exon&& x) = default;
};
inline bool operator==(const Exon& x, const Exon& y) noexcept {
  return std::tie(x.gene, x.coordinates) == std::tie(y.gene, y.coordinates);
}

// override boost hashing
std::size_t hash_value(const Exon& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// specialize std::hash for Exon
namespace std {
template <> struct hash<majiq::Exon> {
  std::size_t operator()(const majiq::Exon& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
/**
 * Exons in sorted contiguous container, with extra set for exons that could
 * not be easily added in sorted order to defer sorting/reordering vector until
 * necessary
 */
template <class CompareExonT = std::less<>>
class Exons : public detail::GeneRegions<Exon, CompareExonT> {
 public:
  using BaseT = detail::GeneRegions<Exon, CompareExonT>;
  using BaseT::remap_genes;
  using BaseT::size;
  using BaseT::make_contiguous;

  Exons() = default;
  Exons(const Exons& x) = default;
  Exons(Exons&& x) = default;
  Exons& operator=(const Exons& x) = default;
  Exons& operator=(Exons&& x) = default;
};
}  // namespace majiq

#endif  // MAJIQ_EXONS_HPP
