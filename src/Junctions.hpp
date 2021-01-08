/**
 * Junctions.hpp
 *
 * Junctions for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_JUNCTIONS_HPP
#define MAJIQ_JUNCTIONS_HPP

#include <tuple>
#include <functional>
#include <boost/functional/hash.hpp>

#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "Connection.hpp"


namespace majiq {
struct ContigJunction : public detail::ContigRegion<OpenInterval> {
 public:
  // constructors
  ContigJunction(KnownContig contig, OpenInterval coordinates,
      GeneStrandness strand)
      : detail::ContigRegion<OpenInterval>{contig, coordinates, strand} {
  }
  ContigJunction(const ContigJunction& x) = default;
  ContigJunction(ContigJunction&& x) = default;
  ContigJunction& operator=(const ContigJunction& x) = default;
  ContigJunction& operator=(ContigJunction&& x) = default;
};
inline bool operator==(
    const ContigJunction& x, const ContigJunction& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    == std::tie(y.contig, y.coordinates, y.strand);
}

struct GeneJunction
  : public detail::GeneRegion<OpenInterval>,
    public detail::Connection {
 public:
  // constructors
  GeneJunction(KnownGene gene, OpenInterval coordinates,
      bool denovo, bool passed_build, bool simplified)
      : detail::GeneRegion<OpenInterval>{gene, coordinates},
        detail::Connection{denovo, passed_build, simplified} {
  }
  GeneJunction(KnownGene gene, OpenInterval coordinates)
      : GeneJunction{gene, coordinates, false, false, false} {
  }
  GeneJunction(const GeneJunction& x) = default;
  GeneJunction(GeneJunction&& x) = default;
  GeneJunction& operator=(const GeneJunction& x) = default;
  GeneJunction& operator=(GeneJunction&& x) = default;
  /**
   * Do we match with specified contig junction (strand matches if same or if
   * contig junction has ambiguous strand)
   */
  // matching with contig junction -- awareness of ambiguous gene strandness
  bool matches(const ContigJunction& rhs) const {
    if (coordinates != rhs.coordinates) {
      return false;
    } else {
      // need to compare information from gene
      const Gene& g = gene.get();
      return g.contig == rhs.contig
        && (rhs.strand == GeneStrandness::AMBIGUOUS || rhs.strand == g.strand);
    }
  }
};
inline bool operator==(const GeneJunction& x, const GeneJunction& y) noexcept {
  return std::tie(x.gene, x.coordinates) == std::tie(y.gene, y.coordinates);
}


// override boost::hash
std::size_t hash_value(const ContigJunction& x) noexcept {
  std::size_t result = hash_value(x.contig);
  boost::hash_combine(result, x.coordinates);
  boost::hash_combine(result, x.strand);
  return result;
}
std::size_t hash_value(const GeneJunction& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// override std::hash for GeneJunction, ContigJunction
namespace std {
template <> struct hash<majiq::ContigJunction> {
  std::size_t operator()(const majiq::ContigJunction& x) const noexcept {
    return majiq::hash_value(x);
  }
};
template <> struct hash<majiq::GeneJunction> {
  std::size_t operator()(const majiq::GeneJunction& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
using GeneJunctions
  = detail::GeneRegions<GeneJunction, std::less<GeneJunction>>;
}  // namespace majiq

#endif  // MAJIQ_JUNCTIONS_HPP
