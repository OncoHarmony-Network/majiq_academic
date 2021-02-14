/**
 * GeneJunctions.hpp
 *
 * Junctions for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_GENEJUNCTIONS_HPP
#define MAJIQ_GENEJUNCTIONS_HPP

#include <tuple>
#include <functional>
#include <boost/functional/hash.hpp>

#include "GeneRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "Connection.hpp"


namespace majiq {

struct GeneJunction
  : public detail::GeneRegion<OpenInterval, detail::Connection> {
 public:
  using BaseT = detail::GeneRegion<OpenInterval, detail::Connection>;
  // access data nicely
  const bool& denovo() const noexcept { return data.denovo; }
  bool& denovo() noexcept { return data.denovo; }
  bool& passed_build() const noexcept { return data.passed_build; }
  bool& passed_build() noexcept { return data.passed_build; }
  bool& simplified() const noexcept { return data.simplified; }
  bool& simplified() noexcept { return data.simplified; }

  // constructors
  GeneJunction(KnownGene gene, OpenInterval coordinates,
      bool denovo, bool passed_build, bool simplified)
      : BaseT{gene, coordinates,
        detail::Connection{denovo, passed_build, simplified}} {
  }
  GeneJunction(KnownGene gene, OpenInterval coordinates)
      : GeneJunction{gene, coordinates, false, false, false} {
  }
  GeneJunction() : GeneJunction{KnownGene{}, OpenInterval{}} { }
  GeneJunction(const GeneJunction& x) = default;
  GeneJunction(GeneJunction&& x) = default;
  GeneJunction& operator=(const GeneJunction& x) = default;
  GeneJunction& operator=(GeneJunction&& x) = default;
};

// override boost::hash
inline std::size_t hash_value(const GeneJunction& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// override std::hash for GeneJunction
namespace std {
template <> struct hash<majiq::GeneJunction> {
  std::size_t operator()(const majiq::GeneJunction& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
using GeneJunctions = detail::Regions<GeneJunction, true>;
}  // namespace majiq

#endif  // MAJIQ_GENEJUNCTIONS_HPP
