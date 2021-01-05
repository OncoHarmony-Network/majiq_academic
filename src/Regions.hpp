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
};

// order regions by genomic position and strand
inline bool operator<(const ContigRegion& x, const ContigRegion& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
/**
 * Order regions by gene, then position
 */
inline bool operator<(const GeneRegion& x, const GeneRegion& y) noexcept {
  return std::tie(x.gene, x.coordinates) < std::tie(y.gene, y.coordinates);
}
/**
 * Compare to genes
 */
inline bool operator<(const GeneRegion& lhs, const KnownGene& rhs) noexcept {
  return lhs.gene < rhs;
}
inline bool operator<(const KnownGene& lhs, const GeneRegion rhs) noexcept {
  return lhs < rhs.gene;
}
/**
 * Compare to contig regions
 */
inline bool operator<(
    const GeneRegion& x, const ContigRegion& y) noexcept {
  const Gene& gx = x.gene.get();
  return std::tie(gx.contig, x.coordinates, gx.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
inline bool operator<(const ContigRegion& x, const GeneRegion& y) noexcept {
  const Gene& gy = y.gene.get();
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(gy.contig, y.coordinates, gy.strand);
}

// unstranded comparison with contig region
inline bool CompareUnstranded(
    const ContigRegion& x, const ContigRegion& y) noexcept {
  return std::tie(x.contig, x.coordinates) < std::tie(y.contig, y.coordinates);
}
inline bool CompareUnstranded(
    const ContigRegion& x, const GeneRegion& y) noexcept {
  const Gene& gy = y.gene.get();
  return std::tie(x.contig, x.coordinates) < std::tie(gy.contig, y.coordinates);
}
inline bool CompareUnstranded(
    const GeneRegion& x, const ContigRegion& y) noexcept {
  const Gene& gx = x.gene.get();
  return std::tie(gx.contig, x.coordinates) < std::tie(y.contig, y.coordinates);
}

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

// derived comparisons (ContigRegion, ContigRegion)
inline bool operator>(const ContigRegion& x, const ContigRegion& y) noexcept {
  return y < x;
}
inline bool operator<=(const ContigRegion& x, const ContigRegion& y) noexcept {
  return !(y < x);
}
inline bool operator>=(const ContigRegion& x, const ContigRegion& y) noexcept {
  return !(x < y);
}

// derived comparisons (GeneRegion, GeneRegion)
inline bool operator>(const GeneRegion& x, const GeneRegion& y) noexcept {
  return y < x;
}
inline bool operator>=(const GeneRegion& x, const GeneRegion& y) noexcept {
  return !(x < y);
}
inline bool operator<=(const GeneRegion& x, const GeneRegion& y) noexcept {
  return !(y < x);
}

// derived comparisons (KnownGene, GeneRegion)
inline bool operator>(const GeneRegion& x, const KnownGene& y) noexcept {
  return y < x;
}
inline bool operator>(const KnownGene& x, const GeneRegion y) noexcept {
  return y < x;
}
inline bool operator<=(const GeneRegion& x, const KnownGene& y) noexcept {
  return !(y < x);
}
inline bool operator<=(const KnownGene& x, const GeneRegion y) noexcept {
  return !(y < x);
}
inline bool operator>=(const GeneRegion& x, const KnownGene& y) noexcept {
  return !(x < y);
}
inline bool operator>=(const KnownGene& x, const GeneRegion y) noexcept {
  return !(x < y);
}

// derived comparisons (ContigRegion, GeneRegion)
inline bool operator>(const GeneRegion& x, const ContigRegion& y) noexcept {
  return y < x;
}
inline bool operator>(const ContigRegion& x, const GeneRegion y) noexcept {
  return y < x;
}
inline bool operator<=(const GeneRegion& x, const ContigRegion& y) noexcept {
  return !(y < x);
}
inline bool operator<=(const ContigRegion& x, const GeneRegion y) noexcept {
  return !(y < x);
}
inline bool operator>=(const GeneRegion& x, const ContigRegion& y) noexcept {
  return !(x < y);
}
inline bool operator>=(const ContigRegion& x, const GeneRegion y) noexcept {
  return !(x < y);
}
}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_REGIONS_HPP
