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
#include <utility>
#include <stdexcept>
#include <tuple>
#include <functional>
#include <algorithm>

#include "Exons.hpp"
#include "GeneJunctions.hpp"
#include "GeneIntrons.hpp"

namespace majiq {

// carry index to either intron or junction
enum class EventType : unsigned char {
  SRC_EVENT = 's',
  DST_EVENT = 't'
};
struct EventIndex {
  size_t ref_exon_idx_;
  EventType event_type_;
};
inline bool operator<(const EventIndex& x, const EventIndex& y) noexcept {
  return std::tie(x.ref_exon_idx_, x.event_type_)
    < std::tie(y.ref_exon_idx_, y.event_type_);
}
inline bool operator==(const EventIndex& x, const EventIndex& y) noexcept {
  return std::tie(x.ref_exon_idx_, x.event_type_)
    == std::tie(y.ref_exon_idx_, y.event_type_);
}
inline bool operator!=(const EventIndex& x, const EventIndex& y) noexcept {
  return !(x == y);
}
struct EventConnection {
  EventIndex event_idx_;
  bool is_intron_;
  size_t connection_idx_;
};
// we later define how to make a comparator of these given junctions/introns

class EventConnections {
 private:
  // pointers to connections
  const std::shared_ptr<GeneJunctions> junctions_;
  const std::shared_ptr<GeneIntrons> introns_;

  // get exons inferred from junctions_ and introns_ (which must not be null)
  const std::shared_ptr<Exons> exons_;
  static std::shared_ptr<Exons> InferredExons(
      const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<GeneIntrons>& introns);
  void check_exons() const {
    if (exons_ != InferredExons(junctions_, introns_)) {
      throw std::runtime_error(
          "EventConnections exons no longer match junctions/introns");
    }
  }

  // track junctions and introns together
  const std::vector<EventConnection> event_connections_;
  static std::vector<EventConnection> MakeSortedEventConnections(
      const GeneJunctions& junctions, const GeneIntrons& introns);

  // track the events themselves and their offsets into event_connections_
  const std::vector<size_t> event_offsets_;
  static std::vector<size_t> CalculateOffsets(
      const std::vector<EventConnection>& event_connections) {
    std::vector<size_t> result = {0};
    auto prev_event_it = event_connections.begin();
    for (auto event_it = event_connections.begin();
        event_it != event_connections.end();
        ++event_it) {
      if (prev_event_it->event_idx_ != event_it->event_idx_) {
        result.push_back(event_it - event_connections.begin());
        prev_event_it = event_it;
      }
    }
    result.push_back(event_connections.size());
    return result;
  }
  // keep the index separately as well
  const std::vector<EventIndex> events_;
  static std::vector<EventIndex> GetIndex(
      const std::vector<EventConnection>& event_connections,
      const std::vector<size_t>& event_offsets) {
    std::vector<EventIndex> result(event_offsets.size() - 1);
    std::transform(event_offsets.begin(), event_offsets.end() - 1,
        result.begin(),
        [&event_connections](size_t i) {
        return event_connections[i].event_idx_; });
    return result;
  }

 public:
  EventConnections(
      const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<GeneIntrons>& introns)
      : junctions_{junctions},
        introns_{introns},
        exons_{InferredExons(junctions_, introns_)},
        event_connections_{MakeSortedEventConnections(*junctions_, *introns_)},
        event_offsets_{CalculateOffsets(event_connections_)},
        events_{GetIndex(event_connections_, event_offsets_)} { }

  const std::vector<EventConnection>& event_connections() const {
    return event_connections_;
  }

  // total number of connections of any kind (including exitrons, simplified)
  size_t size() const noexcept { return event_connections_.size(); }
  // total number of potential events, including redundant and constitutive
  size_t num_potential_events() const noexcept { return events_.size(); }
};

inline std::function<bool(const EventConnection&, const EventConnection&)>
GetEventConnectionComparator(
    const GeneJunctions& junctions, const GeneIntrons& introns) {
  auto result = [&junctions, &introns](
      const EventConnection& x, const EventConnection& y) -> bool {
    if (std::tie(x.event_idx_, x.is_intron_)
        < std::tie(y.event_idx_, y.is_intron_)) {
      return true;
    } else if (std::tie(y.event_idx_, y.is_intron_)
        < std::tie(x.event_idx_, x.is_intron_)) {
      return false;
    } else if (x.connection_idx_ == y.connection_idx_) {
      return false;
    } else {
      // from same event, same type (junction vs intron), and not the same
      auto compare_regions = [](
          // only used when we can assume xr != yr, especially with coordinates
          const auto& xr, const auto& yr, bool is_src) {
        const bool is_forward = xr.gene.strand() == GeneStrandness::FORWARD;
        if (is_src == is_forward) {
          return is_forward == (xr.coordinates < yr.coordinates);
        } else {
          return is_forward == (
              std::tie(xr.coordinates.end, xr.coordinates.start)
              < std::tie(yr.coordinates.end, yr.coordinates.start));
        }
      };
      if (x.is_intron_) {
        throw std::logic_error(
            "EventConnections has reference exon with two distinct introns");
      } else {
        const bool is_src = x.event_idx_.event_type_ == EventType::SRC_EVENT;
        return compare_regions(
            junctions[x.connection_idx_], junctions[y.connection_idx_], is_src);
      }
    }
  };
  return result;
}

inline std::vector<EventConnection>
EventConnections::MakeSortedEventConnections(
    const GeneJunctions& junctions, const GeneIntrons& introns) {
  // initialize results
  std::vector<EventConnection> result;
  result.reserve(2 * (junctions.size() + introns.size()));
  // copy over indexes
  for (size_t i = 0; i < junctions.size(); ++i) {
    const GeneJunction& junction = junctions[i];
    result.push_back(
        EventConnection{
          EventIndex{junction.src_exon_idx(), EventType::SRC_EVENT},
          false, i});
    result.push_back(
        EventConnection{
          EventIndex{junction.dst_exon_idx(), EventType::DST_EVENT},
          false, i});
  }
  for (size_t i = 0; i < introns.size(); ++i) {
    const GeneIntron& intron = introns[i];
    result.push_back(
        EventConnection{
          EventIndex{intron.src_exon_idx(), EventType::SRC_EVENT},
          true, i});
    result.push_back(
        EventConnection{
          EventIndex{intron.dst_exon_idx(), EventType::DST_EVENT},
          true, i});
  }
  // get comparator to sort this
  auto comparator = GetEventConnectionComparator(junctions, introns);
  std::sort(result.begin(), result.end(), comparator);
  return result;
}
inline std::shared_ptr<Exons> EventConnections::InferredExons(
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

}  // namespace majiq


#endif  // MAJIQ_EVENTCONNECTIONS_HPP
