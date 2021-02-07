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
#include <vector>
#include <boost/functional/hash.hpp>

#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "Connection.hpp"


namespace majiq {
struct Intron : public detail::GeneRegion<ClosedInterval, detail::Connection> {
 public:
  using BaseT = detail::GeneRegion<ClosedInterval, detail::Connection>;
  // access data nicely
  const bool& denovo() const noexcept { return data.denovo; }
  bool& denovo() noexcept { return data.denovo; }
  const bool& passed_build() const noexcept { return data.passed_build; }
  bool& passed_build() noexcept { return data.passed_build; }
  const bool& simplified() const noexcept { return data.simplified; }
  bool& simplified() noexcept { return data.simplified; }

  // constructors
  Intron(KnownGene gene, ClosedInterval coordinates,
      bool denovo, bool passed_build, bool simplified)
      : BaseT{gene, coordinates,
        detail::Connection{denovo, passed_build, simplified}} {
  }
  Intron(KnownGene gene, ClosedInterval coordinates)
      : Intron{gene, coordinates, false, false, false} {
  }
  Intron() : Intron{KnownGene{}, ClosedInterval{}} { }
  Intron(const Intron& x) = default;
  Intron(Intron&& x) = default;
  Intron& operator=(const Intron& x) = default;
  Intron& operator=(Intron&& x) = default;
};
inline bool operator==(const Intron& x, const Intron& y) noexcept {
  return std::tie(x.gene, x.coordinates) == std::tie(y.gene, y.coordinates);
}

// override boost hashing
inline std::size_t hash_value(const Intron& x) noexcept {
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
using Introns = detail::Regions<Intron, false>;

}  // namespace majiq

#endif  // MAJIQ_INTRONS_HPP
