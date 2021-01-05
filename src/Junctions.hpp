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
#include <boost/functional/hash.hpp>

#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
struct ContigJunction : public detail::ContigRegion<OpenInterval> {
 public:
  // constructors
  ContigJunction(KnownContig _contig, OpenInterval _coordinates,
      GeneStrandness _strand)
      : detail::ContigRegion<OpenInterval>{_contig, _coordinates, _strand} {
  }
  ContigJunction(const ContigJunction& x)
      : ContigJunction{x.contig, x.coordinates, x.strand} {
  }

  // comparison for equality
  bool operator==(const ContigJunction& rhs) const {
    return std::tie(contig, coordinates, strand)
      == std::tie(rhs.contig, rhs.coordinates, rhs.strand);
  }
};

struct GeneJunction : public detail::GeneRegion<OpenInterval> {
 public:
  // constructors
  GeneJunction(KnownGene _gene, OpenInterval _coordinates)
      : detail::GeneRegion<OpenInterval>{_gene, _coordinates} {
  }
  GeneJunction(const GeneJunction& x) : GeneJunction{x.gene, x.coordinates} {}

  // comparison for equality
  bool operator==(const GeneJunction& rhs) const {
    return std::tie(gene, coordinates)
      == std::tie(rhs.gene, rhs.coordinates);
  }

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
/**
 * GeneJunctions in sorted contiguous container, with extra set for junctions
 * that could not be easily added in sorted order to defer sorting/reordering
 * vector until necessary
 */
class GeneJunctions : public detail::GeneRegions<GeneJunction> {
};
}  // namespace majiq

#endif  // MAJIQ_JUNCTIONS_HPP
