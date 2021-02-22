/**
 * GeneIntrons.hpp
 *
 * Gene introns for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_GENEINTRONS_HPP
#define MAJIQ_GENEINTRONS_HPP

#include <tuple>
#include <functional>
#include <vector>
#include <boost/functional/hash.hpp>

#include "GeneRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "Connection.hpp"


namespace majiq {
struct GeneIntron : public detail::GeneRegion<ClosedInterval, detail::Connection> {
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
  GeneIntron(KnownGene gene, ClosedInterval coordinates,
      bool denovo, bool passed_build, bool simplified)
      : BaseT{gene, coordinates,
        detail::Connection{denovo, passed_build, simplified}} {
  }
  GeneIntron(KnownGene gene, ClosedInterval coordinates)
      : GeneIntron{gene, coordinates, false, false, false} {
  }
  GeneIntron() : GeneIntron{KnownGene{}, ClosedInterval{}} { }
  GeneIntron(const GeneIntron& x) = default;
  GeneIntron(GeneIntron&& x) = default;
  GeneIntron& operator=(const GeneIntron& x) = default;
  GeneIntron& operator=(GeneIntron&& x) = default;
};

// override boost hashing
inline std::size_t hash_value(const GeneIntron& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// specialize std::hash for GeneIntron
namespace std {
template <> struct hash<majiq::GeneIntron> {
  std::size_t operator()(const majiq::GeneIntron& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
using GeneIntrons = detail::Regions<GeneIntron, false>;

}  // namespace majiq

#endif  // MAJIQ_GENEINTRONS_HPP
