/**
 * Junctions.hpp
 *
 * Junctions for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_JUNCTIONS_HPP
#define MAJIQ_JUNCTIONS_HPP

#include <iostream>
#include <tuple>
#include <boost/functional/hash.hpp>

#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
struct ContigJunction {
 public:
  // location
  KnownContig contig;
  ClosedInterval coordinates;
  GeneStrandness strand;

  // constructors
  ContigJunction(KnownContig _contig, ClosedInterval _coordinates,
      GeneStrandness _strand)
      : contig{_contig}, coordinates{_coordinates}, strand{_strand} {
  }
  ContigJunction(const ContigJunction& x)
      : ContigJunction{x.contig, x.coordinates, x.strand} {
  }

  /**
   * Order junctions by genomic position
   */
  bool operator<(const ContigJunction& rhs) const {
    return std::tie(contig, coordinates, strand)
      < std::tie(rhs.contig, rhs.coordinates, rhs.strand);
  }
  bool operator==(const ContigJunction& rhs) const {
    return std::tie(contig, coordinates, strand)
      == std::tie(rhs.contig, rhs.coordinates, rhs.strand);
  }
};

struct GeneJunction {
  // location
  KnownGene gene;
  ClosedInterval coordinates;

  // constructors
  GeneJunction(KnownGene _gene, ClosedInterval _coordinates)
      : gene{_gene}, coordinates{_coordinates} {
  }
  GeneJunction(const GeneJunction& x) : GeneJunction{x.gene, x.coordinates} {}

  // downcast to junction
  ContigJunction AsContigJunction() const {
    const Gene& g = gene.get();
    return ContigJunction{g.contig, coordinates, g.strand};
  }

  /**
   * Order junctions by genomic position
   */
  bool operator<(const GeneJunction& rhs) const {
    return std::tie(gene, coordinates)
      < std::tie(rhs.gene, rhs.coordinates);
  }
  bool operator==(const GeneJunction& rhs) const {
    return std::tie(gene, coordinates)
      == std::tie(rhs.gene, rhs.coordinates);
  }
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

// allow Junction to be passed into output stream (e.g. std::cout)
std::ostream& operator<<(std::ostream& os, const ContigJunction& x) noexcept {
  os << x.contig.get()
    << ":" << x.strand
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
}
std::ostream& operator<<(std::ostream& os, const GeneJunction& x) noexcept {
  os << x.gene.get()
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
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
#endif  // MAJIQ_JUNCTIONS_HPP
