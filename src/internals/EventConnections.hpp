/**
 * EventConnections.hpp
 *
 * Identify and appropriately sort connections by src/dst (stranded) exons for
 * the purposes of simplifying and identifying events (splicing and
 * constitutive)
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_EVENTCONNECTIONS_HPP
#define MAJIQ_EVENTCONNECTIONS_HPP

#include <memory>
#include <vector>
#include <stdexcept>
#include <algorithm>

#include "Exons.hpp"
#include "GeneJunctions.hpp"
#include "GeneIntrons.hpp"
#include "EventConnection.hpp"

namespace majiq {

class EventConnections {
 private:
  // pointers to connections
  const std::shared_ptr<GeneJunctions> junctions_;
  const std::shared_ptr<GeneIntrons> introns_;
  const std::shared_ptr<Exons> exons_;  // needs to match junctions and introns
  const std::vector<EventConnection> event_connections_;

 public:
  size_t size() const noexcept { return event_connections_.size(); }
  const std::shared_ptr<Exons>& exons() const { return exons_; }
  const std::shared_ptr<GeneJunctions>& junctions() const { return junctions_; }
  const std::shared_ptr<GeneIntrons>& introns() const { return introns_; }
  const std::vector<EventConnection> event_connections() const {
    return event_connections_;
  }

  EventConnectionReference operator[](size_t ec_idx) const {
    return EventConnectionReference{
        event_connections_[ec_idx], *exons_, *junctions_, *introns_};
  }

 private:
  static std::shared_ptr<Exons> InferredExons(
      const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<GeneIntrons>& introns) {
    if (junctions == nullptr) {
      throw std::runtime_error(
          "EventConnections junctions cannot be null-valued");
    } else if (introns == nullptr) {
      throw std::runtime_error(
          "EventConnections introns cannot be null-valued");
    } else if (!(
          junctions->is_connected()
          && introns->is_connected()
          && junctions->connected_exons() == introns->connected_exons())) {
      throw std::runtime_error(
          "EventConnections junctions/introns not connected to same exons");
    }
    return junctions->connected_exons();
  }

  static std::vector<EventConnection> MakeSortedEventConnections(
      const GeneJunctions& junctions, const GeneIntrons& introns) {
    // initialize results
    std::vector<EventConnection> result;
    result.reserve(2 * (junctions.size() + introns.size()));
    // copy over indexes
    for (size_t i = 0; i < junctions.size(); ++i) {
      const GeneJunction& junction = junctions[i];
      result.push_back(
          EventConnection{
            Event{junction.src_exon_idx(), EventType::SRC_EVENT},
            false, i});
      result.push_back(
          EventConnection{
            Event{junction.dst_exon_idx(), EventType::DST_EVENT},
            false, i});
    }
    for (size_t i = 0; i < introns.size(); ++i) {
      const GeneIntron& intron = introns[i];
      result.push_back(
          EventConnection{
            Event{intron.src_exon_idx(), EventType::SRC_EVENT},
            true, i});
      result.push_back(
          EventConnection{
            Event{intron.dst_exon_idx(), EventType::DST_EVENT},
            true, i});
    }
    std::sort(result.begin(), result.end(),
        [&junctions, &introns, &exons = *junctions.connected_exons()](
          const EventConnection& x, const EventConnection& y) -> bool {
        return (
            EventConnectionReference{x, exons, junctions, introns}
            < EventConnectionReference{y, exons, junctions, introns}); });
    return result;
  }

  // since junctions, introns mutable with connections, check if not match
  void check_exons() const {
    if (exons_ != InferredExons(junctions_, introns_)) {
      throw std::runtime_error(
          "EventConnections exons no longer match junctions/introns");
    }
  }

 public:
  EventConnections(
      const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<GeneIntrons>& introns)
      : junctions_{junctions},
        introns_{introns},
        exons_{InferredExons(junctions_, introns_)},
        event_connections_{
          MakeSortedEventConnections(*junctions_, *introns_)} { }
};

}  // namespace majiq


#endif  // MAJIQ_EVENTCONNECTIONS_HPP
