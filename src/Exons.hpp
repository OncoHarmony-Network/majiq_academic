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
  using BaseT = detail::GeneRegion<ClosedInterval>;
  struct DefaultAnnotated { };
  struct MakeDenovo { };

  // exon-specific members
  ClosedInterval annotated_coordinates;

  // exon-specific information
  inline bool is_denovo() const noexcept {
    return annotated_coordinates.is_invalid();
  }
  inline bool is_half_exon() const noexcept {
    return coordinates.is_half_interval();
  }
  /**
   * If denovo, return current coordinates to match previous MAJIQ
   */
  inline const ClosedInterval& legacy_annotated_coordinates() const noexcept {
    return is_denovo() ? coordinates : annotated_coordinates;
  }

  // constructors
  Exon(KnownGene _gene, ClosedInterval _coordinates, ClosedInterval _annotated)
      : BaseT{_gene, _coordinates}, annotated_coordinates{_annotated} {
  }
  Exon(KnownGene _gene, ClosedInterval _coordinates, DefaultAnnotated)
      : Exon{_gene, _coordinates, _coordinates} {
  }
  Exon(KnownGene _gene, ClosedInterval _coordinates, MakeDenovo)
      : Exon{_gene, _coordinates, ClosedInterval{-1, -1}} {
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
using Exons = detail::GeneRegions<Exon, std::less<Exon>>;
}  // namespace majiq

#endif  // MAJIQ_EXONS_HPP
