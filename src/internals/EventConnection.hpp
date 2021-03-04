/**
 * EventConnection.hpp
 *
 * Representation of events and connections that are part of them
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_EVENTCONNECTION_HPP
#define MAJIQ_EVENTCONNECTION_HPP

#include <tuple>
#include <variant>

#include "MajiqTypes.hpp"
#include "Exons.hpp"
#include "GeneIntrons.hpp"
#include "GeneJunctions.hpp"


namespace majiq {

struct EventConnection {
  Event event_;
  bool is_intron_;
  size_t connection_idx_;
};

class EventConnectionReference {
 private:
  const EventConnection& raw_;
  const Exons& exons_;
  const GeneJunctions& junctions_;
  const GeneIntrons& introns_;

 public:
  EventConnectionReference(
      const EventConnection& raw,
      const Exons& exons,
      const GeneJunctions& junctions,
      const GeneIntrons& introns)
      : raw_{raw},
        exons_{exons},
        junctions_{junctions},
        introns_{introns} { }

  const EventConnection& raw() const { return raw_; }
  const Event& event() const { return raw().event_; }
  const bool& is_intron() const { return raw().is_intron_; }
  const bool is_junction() const { return !is_intron(); }
  const size_t& connection_idx() const { return raw().connection_idx_; }

  const size_t& ref_exon_idx() const { return event().ref_exon_idx_; }
  const EventType& type() const { return event().type_; }
  bool is_source() const { return type() == EventType::SRC_EVENT; }
  const Exon& ref_exon() const { return exons_[ref_exon_idx()]; }
  const KnownGene& gene() const { return ref_exon().gene; }
  GeneStrandness strand() const { return gene().strand(); }
  bool is_forward() const { return strand() == GeneStrandness::FORWARD; }

  bool ref_is_start() const { return is_source() == is_forward(); }

  // If we know for sure that it's a junction or intron
  const GeneJunction& junction() const {
    return junctions_[connection_idx()];
  }
  const GeneIntron& intron() const { return introns_[connection_idx()]; }

  // otherwise, we need to handle both cases
  using connection_t = std::variant<GeneJunction, GeneIntron>;
  connection_t connection() const {
    return is_intron()
      ? connection_t{intron()} : connection_t{junction()};
  }

  const position_t& start() const {
    return is_intron()
      ? intron().coordinates.start : junction().coordinates.start;
  }
  const position_t& end() const {
    return is_intron()
      ? intron().coordinates.end : junction().coordinates.end;
  }
  const position_t& ref_coordinate() const {
    return ref_is_start() ? start() : end();
  }
  const position_t& other_coordinate() const {
    return ref_is_start() ? end() : start();
  }

  const size_t& other_exon_idx() const {
    return ref_is_start()
      ? (is_intron() ? intron().end_exon_idx() : junction().end_exon_idx())
      : (is_intron() ? intron().start_exon_idx() : junction().start_exon_idx());
  }
  const Exon& other_exon() const { return exons_[other_exon_idx()]; }

  const bool& simplified() const {
    return is_intron() ? intron().simplified() : junction().simplified();
  }
  bool passed_build() const {
    return is_intron() ? intron().passed_build() : junction().passed_build();
  }
  bool is_exitron() const {
    return is_intron() ? intron().is_exitron() : junction().is_exitron();
  }
  bool for_event() const { return !(is_exitron() || simplified()); }
  bool for_passed() const { return passed_build() && for_event(); }
};

inline bool operator<(
    const EventConnectionReference& x, const EventConnectionReference& y) {
  // easy if different events or if junction vs intron
  if (std::tie(x.event(), x.is_intron())
      != std::tie(y.event(), y.is_intron())) {
    return std::tie(x.event(), x.is_intron())
      < std::tie(y.event(), y.is_intron());
  }
  // otherwise, same event, and both are junction or intron
  if (x.connection_idx() == y.connection_idx()) {
    // rule out comparing to same event connection
    return false;
  } else if (x.is_intron()) {
    // x and y are from same event, are both introns, but not same connection
    throw std::logic_error(
        "EventConnections has event with two distinct introns");
  } else {
    // comparing junctions
    const OpenInterval& xc = x.junction().coordinates;
    const OpenInterval& yc = y.junction().coordinates;
    const bool is_source = x.is_source();
    const bool is_forward = x.is_forward();
    if (is_source == is_forward) {
      return is_forward == (xc < yc);
    } else {
      return is_forward
        == (std::tie(xc.end, xc.start) < std::tie(yc.end, yc.start));
    }
  }
}


}  // namespace majiq


#endif  // MAJIQ_EVENTCONNECTION_HPP
