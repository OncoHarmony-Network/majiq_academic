/**
 * Introns.hpp
 *
 * Introns for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_INTRONS_HPP
#define MAJIQ_INTRONS_HPP

#include <tuple>
#include <functional>
#include <boost/functional/hash.hpp>

#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
struct Intron : detail::GeneRegion<ClosedInterval> {
 public:
  // constructors
  Intron(KnownGene _gene, ClosedInterval _coordinates)
      : detail::GeneRegion<ClosedInterval>{_gene, _coordinates} {
  }
  Intron(const Intron& x) = default;
  Intron(Intron&& x) = default;
  Intron& operator=(const Intron& x) = default;
  Intron& operator=(Intron&& x) = default;
};
inline bool operator==(const Intron& x, const Intron& y) noexcept {
  return std::tie(x.gene, x.coordinates) == std::tie(y.gene, y.coordinates);
}

// override boost hashing
std::size_t hash_value(const Intron& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// specialize std::hash for Intron
namespace std {
template <> struct hash<majiq::Intron> {
  std::size_t operator()(const majiq::Intron& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
using Introns = detail::GeneRegions<Intron, std::less<Intron>>;
}  // namespace majiq

#endif  // MAJIQ_INTRONS_HPP
