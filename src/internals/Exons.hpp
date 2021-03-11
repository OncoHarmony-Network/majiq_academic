/**
 * Exons.hpp
 *
 * Exons for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_EXONS_HPP
#define MAJIQ_EXONS_HPP

#include <tuple>
#include <functional>
#include <vector>
#include <algorithm>
#include <utility>
#include <memory>
#include <boost/functional/hash.hpp>

#include "GeneRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "checksum.hpp"


namespace majiq {
struct Exon : detail::GeneRegion<ClosedInterval, ClosedInterval> {
 public:
  using BaseT = detail::GeneRegion<ClosedInterval, ClosedInterval>;
  struct DefaultAnnotated { };
  struct MakeDenovo { };

  // exon-specific members
  inline ClosedInterval& annotated_coordinates() noexcept { return data; }
  inline const ClosedInterval& annotated_coordinates() const noexcept {
    return data;
  }

  // constructors
  Exon(KnownGene _gene, ClosedInterval _coordinates, ClosedInterval _annotated)
      : BaseT{_gene, _coordinates, _annotated} {
  }
  Exon(KnownGene _gene, ClosedInterval _coordinates, DefaultAnnotated)
      : Exon{_gene, _coordinates, _coordinates} {
  }
  Exon(KnownGene _gene, ClosedInterval _coordinates, MakeDenovo)
      : Exon{_gene, _coordinates, ClosedInterval{}} {
  }
  // if no specifier passed, annotated by default
  Exon(KnownGene _gene, ClosedInterval _coordinates)
      : Exon{_gene, _coordinates, DefaultAnnotated{}} { }
  Exon() : Exon{KnownGene{}, ClosedInterval{}, ClosedInterval{}} { }
  Exon(const Exon& x) = default;
  Exon(Exon&& x) = default;
  Exon& operator=(const Exon& x) = default;
  Exon& operator=(Exon&& x) = default;

  // exon-specific information
  inline bool is_denovo() const noexcept {
    return annotated_coordinates().is_invalid();
  }
  inline bool is_full_exon() const noexcept {
    return coordinates.is_full_interval();
  }
  inline bool is_half_exon() const noexcept {
    return coordinates.is_half_interval();
  }
  inline bool is_invalid() const noexcept {
    return coordinates.is_invalid();
  }
  /**
   * If denovo, return current coordinates to match previous MAJIQ
   */
  inline const ClosedInterval& legacy_annotated_coordinates() const noexcept {
    return is_denovo() ? coordinates : annotated_coordinates();
  }
  /**
   * Get annotated version of exon
   *
   * @note if exon is not annotated -- returns exon with same coordinates which
   * is marked as annotated with same coordinates
   */
  Exon get_annotated() const {
    return Exon{gene, legacy_annotated_coordinates(), DefaultAnnotated{}};
  }
};

// override boost hashing
inline std::size_t hash_value(const Exon& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// specialize std::hash for Exon
namespace std {
template <> struct hash<majiq::Exon> {
  std::size_t operator()(const majiq::Exon& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {

class Exons : public detail::Regions<Exon, false> {
  using BaseT = detail::Regions<Exon, false>;

 public:
  Exons(
      const std::shared_ptr<Genes>& genes, std::vector<Exon>&& x)
      : BaseT{genes, std::move(x)} {
    if (parents() == nullptr) {
      throw std::invalid_argument("Exons cannot have null genes");
    }
  }
};

inline detail::checksum_t checksum(const Exons& x) noexcept {
  detail::checksum_gen_t gen;
  for (const auto& e : x) {
    gen.process_bytes(&e.gene.idx_, sizeof(e.gene.idx_));
    gen.process_bytes(&e.coordinates.start, sizeof(e.coordinates.start));
    gen.process_bytes(&e.coordinates.end, sizeof(e.coordinates.end));
    gen.process_bytes(&e.data.start, sizeof(e.data.start));
    gen.process_bytes(&e.data.end, sizeof(e.data.end));
  }
  return detail::checksum_t{gen.checksum()};
}
}  // namespace majiq

#endif  // MAJIQ_EXONS_HPP
