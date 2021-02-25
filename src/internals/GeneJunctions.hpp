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
#include <vector>
#include <memory>
#include <functional>
#include <utility>
#include <boost/functional/hash.hpp>

#include "GeneRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "Connection.hpp"
#include "Exons.hpp"


namespace majiq {

struct GeneJunction
    : public detail::GeneRegion<OpenInterval, detail::Connection> {
 public:
  using BaseT = detail::GeneRegion<OpenInterval, detail::Connection>;
  // access data nicely
  const bool& denovo() const noexcept { return data.denovo; }
  bool& denovo() noexcept { return data.denovo; }
  bool& passed_build() const noexcept { return data.passed_build; }
  bool& simplified() const noexcept { return data.simplified; }
  size_t& start_exon_idx() const noexcept { return data.start_exon_idx; }
  size_t& end_exon_idx() const noexcept { return data.end_exon_idx; }

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
class GeneJunctions : public detail::Regions<GeneJunction, true> {
  using BaseT = detail::Regions<GeneJunction, true>;
 public:
  explicit GeneJunctions(
      const std::shared_ptr<Genes>& genes, std::vector<GeneJunction>&& x)
      : BaseT{genes, std::move(x)} { }
  GeneJunctions() : BaseT{} { }
  GeneJunctions(const GeneJunctions&) = default;
  GeneJunctions(GeneJunctions&&) = default;
  GeneJunctions& operator=(const GeneJunctions&) = default;
  GeneJunctions& operator=(GeneJunctions&&) = default;
};
}  // namespace majiq

#endif  // MAJIQ_GENEJUNCTIONS_HPP
