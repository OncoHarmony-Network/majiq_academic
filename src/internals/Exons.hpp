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
#include <sstream>
#include <string>

#include "GeneRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "checksum.hpp"


namespace majiq {

using ExonIntervalT = ClosedOrHalfInterval;

struct Exon : detail::GeneRegion<ExonIntervalT, ExonIntervalT> {
 public:
  using BaseT = detail::GeneRegion<ExonIntervalT, ExonIntervalT>;
  struct DefaultAnnotated { };
  struct MakeDenovo { };

  // exon-specific members
  inline ExonIntervalT& annotated_coordinates() noexcept { return data; }
  inline const ExonIntervalT& annotated_coordinates() const noexcept {
    return data;
  }
  // event id for this exon with given event type
  std::string event_id(const EventType event_type) const {
    std::ostringstream oss;
    oss << gene.get() << ':' << event_type << ':' << coordinates;
    return oss.str();
  }

  // constructors
  Exon(KnownGene _gene, ExonIntervalT _coordinates, ExonIntervalT _annotated)
      : BaseT{_gene, _coordinates, _annotated} {
  }
  Exon(KnownGene _gene, ExonIntervalT _coordinates, DefaultAnnotated)
      : Exon{_gene, _coordinates, _coordinates} {
  }
  Exon(KnownGene _gene, ExonIntervalT _coordinates, MakeDenovo)
      : Exon{_gene, _coordinates, ExonIntervalT{}} {
  }
  // if no specifier passed, annotated by default
  Exon(KnownGene _gene, ExonIntervalT _coordinates)
      : Exon{_gene, _coordinates, DefaultAnnotated{}} { }
  Exon() : Exon{KnownGene{}, ExonIntervalT{}, ExonIntervalT{}} { }
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
  inline bool is_exon_extension() const noexcept {
    return !is_denovo() && (coordinates != annotated_coordinates());
  }
  /**
   * If denovo, return current coordinates to match previous MAJIQ
   */
  inline const ExonIntervalT& legacy_annotated_coordinates() const noexcept {
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

  /**
   * Return original set of annotated exons (remove denovos, reset coordinates)
   */
  Exons get_annotated() const {
    std::vector<Exon> exon_vec;
    for (auto ex_it = begin(); ex_it != end(); ++ex_it) {
      if (!ex_it->is_denovo()) { exon_vec.push_back(ex_it->get_annotated()); }
    }
    return Exons{parents(), std::move(exon_vec)};
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
