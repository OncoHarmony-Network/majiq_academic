/**
 * Regions.hpp
 *
 * Regions base class for splicegraph (closed intervals, relative to known
 * contigs or known genes)
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_REGIONS_HPP
#define MAJIQ_REGIONS_HPP

#include <iostream>
#include <tuple>

#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
namespace detail {

struct ContigRegion {
 public:
  // location
  KnownContig contig;
  ClosedInterval coordinates;
  GeneStrandness strand;

  // constructors
  ContigRegion(KnownContig _contig, ClosedInterval _coordinates,
      GeneStrandness _strand)
      : contig{_contig}, coordinates{_coordinates}, strand{_strand} {
  }
  ContigRegion(const ContigRegion& x)
      : ContigRegion{x.contig, x.coordinates, x.strand} {
  }

};
// order regions by genomic position
inline bool operator<(const ContigRegion& x, const ContigRegion& y) noexcept {
  return std::tie(x.contig, x.coordinates) < std::tie(y.contig, y.coordinates);
}

struct GeneRegion {
 public:
  // location
  KnownGene gene;
  ClosedInterval coordinates;

  // constructors
  GeneRegion(KnownGene _gene, ClosedInterval _coordinates)
      : gene{_gene}, coordinates{_coordinates} {
  }
  GeneRegion(const GeneRegion& x) : GeneRegion{x.gene, x.coordinates} {}

  /**
   * Order regions by gene, then position
   */
  friend bool operator<(
      const GeneRegion& lhs, const GeneRegion& rhs) noexcept {
    return std::tie(lhs.gene, lhs.coordinates)
      < std::tie(rhs.gene, rhs.coordinates);
  }
  /**
   * Compare to contig regions (ignore strand)
   */
  friend bool operator<(
      const GeneRegion& lhs, const ContigRegion& rhs) noexcept {
    const Gene& g_lhs = lhs.gene.get();
    return std::tie(g_lhs.contig, lhs.coordinates)
      < std::tie(rhs.contig, rhs.coordinates);
  }
  friend bool operator<(
      const ContigRegion& lhs, const GeneRegion& rhs) noexcept {
    const Gene& g_rhs = rhs.gene.get();
    return std::tie(lhs.contig, lhs.coordinates)
      < std::tie(g_rhs.contig, rhs.coordinates);
  }
};

// allow regions to be passed into output stream (e.g. std::cout)
std::ostream& operator<<(std::ostream& os, const ContigRegion& x) noexcept {
  os << x.contig.get()
    << ":" << x.strand
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
}
std::ostream& operator<<(std::ostream& os, const GeneRegion& x) noexcept {
  os << x.gene.get()
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
}


}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_REGIONS_HPP
