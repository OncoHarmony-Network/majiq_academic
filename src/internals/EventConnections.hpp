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
#include <variant>
#include <set>

#include "Exons.hpp"
#include "GeneJunctions.hpp"
#include "GeneIntrons.hpp"

namespace majiq {

// carry index to either intron or junction
enum class EventType : unsigned char {
  SRC_EVENT = 's',
  DST_EVENT = 't'
};
inline EventType OtherEventType(EventType x) {
  switch (x) {
    case EventType::SRC_EVENT:
      return EventType::DST_EVENT;
    case EventType::DST_EVENT:
      return EventType::SRC_EVENT;
    default:
      return x;
  }
}
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

 public:
  EventConnections(
      const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<GeneIntrons>& introns)
      : junctions_{junctions},
        introns_{introns},
        exons_{InferredExons(junctions_, introns_)},
        event_connections_{
          MakeSortedEventConnections(*junctions_, *introns_)} { }

  const std::vector<EventConnection>& event_connections() const {
    return event_connections_;
  }
  const EventConnection& operator[](size_t ec_idx) const {
    return event_connections_[ec_idx];
  }
  std::variant<GeneJunction, GeneIntron> junction_or_intron(
      const EventConnection& ec) const {
    using result_t = std::variant<GeneJunction, GeneIntron>;
    return ec.is_intron_
      ? result_t{(*introns_)[ec.connection_idx_]}
      : result_t{(*junctions_)[ec.connection_idx_]};
  }
  std::variant<GeneJunction, GeneIntron> junction_or_intron(
      size_t ec_idx) const {
    return junction_or_intron((*this)[ec_idx]);
  }
  size_t other_exon_idx(size_t ec_idx) const {
    const EventConnection& ec = (*this)[ec_idx];
    return std::visit(
        [ref_src = ec.event_idx_.event_type_ == EventType::SRC_EVENT](
          const auto& x) -> size_t {
        return ref_src ? x.dst_exon_idx() : x.src_exon_idx(); },
        junction_or_intron(ec));
  }
  bool simplified(size_t ec_idx) const {
    return std::visit(
        [](const auto& x) -> bool { return x.simplified(); },
        junction_or_intron(ec_idx));
  }
  bool passed_build(size_t ec_idx) const {
    return std::visit(
        [](const auto& x) -> bool { return x.passed_build(); },
        junction_or_intron(ec_idx));
  }
  bool is_exitron(size_t ec_idx) const {
    return std::visit(
        [](const auto& x) -> bool { return x.is_exitron(); },
        junction_or_intron(ec_idx));
  }
  // for event: not simplified or exitron
  bool for_event(size_t ec_idx) const {
    return std::visit(
        [](const auto& x) -> bool {
        return !(x.is_exitron() || x.simplified()); },
        junction_or_intron(ec_idx));
  }
  // for passed: not simplified or exitron, but also passed
  bool for_passed(size_t ec_idx) const {
    return std::visit(
        [](const auto& x) -> bool {
        return x.passed_build() && !(x.is_exitron() || x.simplified()); },
        junction_or_intron(ec_idx));
  }
  KnownGene gene(size_t ec_idx) const {
    return std::visit(
        [](const auto& x) -> KnownGene { return x.gene; },
        junction_or_intron(ec_idx));
  }
  position_t coordinates_start(size_t ec_idx) const {
    return std::visit(
        [](const auto& x) -> position_t { return x.coordinates.start; },
        junction_or_intron(ec_idx));
  }
  position_t coordinates_end(size_t ec_idx) const {
    return std::visit(
        [](const auto& x) -> position_t { return x.coordinates.end; },
        junction_or_intron(ec_idx));
  }


  // total number of connections of any kind (including exitrons, simplified)
  size_t size() const noexcept { return event_connections_.size(); }
};

class Events {
 private:
  const std::shared_ptr<EventConnections> event_connections_;
  const std::vector<size_t> offsets_;
  const std::vector<EventIndex> events_;

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
  Events(const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<GeneIntrons>& introns)
      : event_connections_{
          std::make_shared<EventConnections>(junctions, introns)},
        offsets_{CalculateOffsets(event_connections_->event_connections())},
        events_{GetIndex(event_connections_->event_connections(), offsets_)} { }

  const std::shared_ptr<EventConnections>& event_connections() const {
    return event_connections_;
  }
  const std::vector<size_t>& offsets() const { return offsets_; }
  const std::vector<EventIndex>& events() const { return events_; }
  const EventIndex& operator[](size_t event_idx) const {
    return events_[event_idx];
  }
  size_t size() const noexcept { return events_.size(); }
  KnownGene gene(size_t event_idx) const {
    return event_connections_->gene(offsets_[event_idx]);
  }
  bool has_intron(size_t event_idx) const {
    const auto v_connection
      = event_connections_->junction_or_intron(offsets_[1 + event_idx] - 1);
    const auto intron_ptr = std::get_if<GeneIntron>(&v_connection);
    return !(intron_ptr == nullptr || intron_ptr->simplified());
  }
  size_t event_size(size_t event_idx) const {
    size_t ct = 0;
    for (size_t ec_idx = offsets_[event_idx];
        ec_idx < offsets_[1 + event_idx];
        ++ec_idx) {
      if (event_connections_->for_event(ec_idx)) { ++ct; }
    }
    return ct;
  }
  bool passed(size_t event_idx) const {
    for (size_t ec_idx = offsets_[event_idx];
        ec_idx < offsets_[1 + event_idx];
        ++ec_idx) {
      if (event_connections_->for_passed(ec_idx)) { return true; }
    }
    return false;
  }
  std::set<size_t> other_exon_idx_set(size_t event_idx) const {
    std::set<size_t> result;
    for (size_t ec_idx = offsets_[event_idx];
        ec_idx < offsets_[1 + event_idx];
        ++ec_idx) {
      if (event_connections_->for_event(ec_idx)) {
        result.insert(event_connections_->other_exon_idx(ec_idx));
      }
    }
    return result;
  }
  bool redundant(size_t event_idx) const {
    std::set<size_t> other = other_exon_idx_set(event_idx);
    if (other.size() > 1) {
      return false;
    } else if (events_[event_idx].event_type_ == EventType::DST_EVENT) {
      return true;
    } else {
      // only have one exon, but it's a source exon.
      // If the corresponding exon only has one source (the current exon), then
      // they are mutually redundant, and we consider this exon to be not
      // redundant.
      const size_t other_exon_idx = *other.begin();
      const EventIndex other_event = {other_exon_idx, EventType::DST_EVENT};
      // we need to find this other event, which we can do by searching locally
      // we can bound this other event to basically 2 * difference in exon idx
      const size_t bound = 2 + 2 * (
          events_[event_idx].ref_exon_idx_ > other_exon_idx
          ?  events_[event_idx].ref_exon_idx_ - other_exon_idx
          : other_exon_idx - events_[event_idx].ref_exon_idx_);
      const size_t other_lb = event_idx > bound ? event_idx - bound : 0;
      const size_t other_ub = std::min(event_idx + bound, events_.size());
      // we expect it at:
      auto other_it = std::lower_bound(
          events_.begin() + other_lb, events_.begin() + other_ub, other_event);
      // if can't find other event or if it's redundant, then current is not
      return !(other_it == events_.end()
          || *other_it != other_event
          || redundant(other_it - events_.begin()));
    }
  }

  // TODO(jaicher): lsv_id (need to retrieve reference exon)
  // TODO(jaicher): lsv_description
  // TODO(jaicher): simplifier
  // TODO(jaicher): identify true events, LSVs, constitutive
  // TODO(jaicher): filter to LSVs, etc.
  // TODO(jaicher): assign coverage to events
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
