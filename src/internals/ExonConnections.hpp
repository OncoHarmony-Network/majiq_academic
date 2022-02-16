/**
 * ExonConnections.hpp
 *
 * Efficient lookup of junctions and introns associated with each exon in a
 * direction
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_EXONCONNECTIONS_HPP
#define MAJIQ_EXONCONNECTIONS_HPP

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "MajiqConstants.hpp"
#include "MajiqTypes.hpp"
#include "Exons.hpp"
#include "GeneJunctions.hpp"
#include "GeneIntrons.hpp"
#include "Events.hpp"

namespace majiq {

namespace detail {
struct ExonConnectionsIndexes {
  const std::vector<size_t> idx_;
  const std::vector<size_t> exon_offsets_;

  size_t num_connections() const { return idx_.size(); }
  size_t num_exons() const { return exon_offsets_.size() - 1; }

  typename std::vector<size_t>::const_iterator
  begin() const { return idx_.cbegin(); }
  typename std::vector<size_t>::const_iterator
  end() const { return idx_.cend(); }

  typename std::vector<size_t>::const_iterator
  begin_exon(size_t idx) const { return idx_.cbegin() + exon_offsets_[idx]; }
  typename std::vector<size_t>::const_iterator
  end_exon(size_t idx) const { return idx_.cbegin() + exon_offsets_[1 + idx]; }

  ExonConnectionsIndexes(
      std::vector<size_t>&& idx, std::vector<size_t>&& exon_offsets)
      : idx_{std::move(idx)}, exon_offsets_{std::move(exon_offsets)} {
    if (exon_offsets_.empty()) {
      throw std::invalid_argument(
          "Empty offsets in in ExonConnectionsIndexes");
    } else if (idx_.size() != exon_offsets_.back()) {
      throw std::invalid_argument(
          "Offsets do not match connections in ExonConnectionsIndexes");
    }
  }

  template <typename ConnectionsT>
  static ExonConnectionsIndexes FromConnectionsAndExons(
      const ConnectionsT& connections, const std::shared_ptr<Exons>& exons_ptr,
      EventType type) {
    if (!connections.is_connected()) {
      std::runtime_error("Connections unconnected in ExonConnectionsIndexes");
    } else if (connections.connected_exons() != exons_ptr) {
      std::runtime_error(
          "Connections connected to different exons in ExonConnectionsIndexes");
    }
    const Exons& exons = *exons_ptr;
    // construct sorted indexes to connections
    std::vector<size_t> idx(connections.size());
    std::vector<size_t> offsets(1 + exons.size());
    size_t offset_exon_idx = 0;  // current exon for offsets
    for (const KnownGene& gene : *(connections.parents())) {
      const bool strand_forward = gene.strand() == GeneStrandness::FORWARD;
      const bool ref_is_start
        = strand_forward == (type == EventType::SRC_EVENT);
      auto get_ref_exon_idx
        = [&ref_is_start](const auto& x) {
          return ref_is_start ? x.start_exon_idx() : x.end_exon_idx(); };
      auto coordinate_subpredicate
        = [&ref_is_start, &strand_forward](const auto& x, const auto& y) {
          auto xt = ref_is_start
            ? x.coordinates.as_tuple() : x.coordinates.rev_tuple();
          auto yt = ref_is_start
            ? y.coordinates.as_tuple() : y.coordinates.rev_tuple();
          return strand_forward ? xt < yt : yt < xt;
        };
      auto sort_predicate
        = [&connections, &get_ref_exon_idx, &coordinate_subpredicate](
            const size_t& x_idx, const size_t& y_idx) {
          const auto& x = connections[x_idx];
          const auto& y = connections[y_idx];
          size_t xref = get_ref_exon_idx(x);
          size_t yref = get_ref_exon_idx(y);
          return xref < yref || (xref == yref && coordinate_subpredicate(x, y));
        };
      // first/last for the gene
      const auto idx_first
        = idx.begin() + (connections.begin_parent(gene) - connections.begin());
      const auto idx_last
        = idx.begin() + (connections.end_parent(gene) - connections.begin());
      // set up indexes, then sort them for this gene
      std::iota(idx_first, idx_last, idx_first - idx.begin());
      std::sort(idx_first, idx_last, sort_predicate);
      // now get exon offsets
      for (auto it = idx_first; it != idx_last; ++it) {
        const size_t it_ref_exon_idx = get_ref_exon_idx(connections[*it]);
        for (; offset_exon_idx < it_ref_exon_idx; ++offset_exon_idx) {
          offsets[1 + offset_exon_idx] = it - idx.begin();
        }  // update offset_exon_idx to current exon, set value appropriately
      }  // loop over the sorted connection indexes
    }  // done iterating over all genes
    for (; offset_exon_idx < exons.size(); ++offset_exon_idx) {
      offsets[1 + offset_exon_idx] = connections.size();
    }
    return ExonConnectionsIndexes{std::move(idx), std::move(offsets)};
  }
};
}  // namespace detail

class ExonConnections {
 private:
  const std::shared_ptr<Exons> exons_;
  const std::shared_ptr<GeneIntrons> introns_;
  const std::shared_ptr<GeneJunctions> junctions_;

  const detail::ExonConnectionsIndexes src_junctions_;
  const detail::ExonConnectionsIndexes dst_junctions_;
  const detail::ExonConnectionsIndexes src_introns_;
  const detail::ExonConnectionsIndexes dst_introns_;

  const detail::ExonConnectionsIndexes& junctions_for(
      const Event& event) const {
    return event.type_ == EventType::SRC_EVENT
      ? src_junctions_ : dst_junctions_;
  }
  const detail::ExonConnectionsIndexes& introns_for(const Event& event) const {
    return event.type_ == EventType::SRC_EVENT
      ? src_introns_ : dst_introns_;
  }

 public:
  size_t num_exons() const { return exons_->size(); }

  // indicate if event direction is valid. This rules out half exons in
  // directions of missing coordinates
  bool valid_event(const Event& event) const {
    return (*exons_)[event.ref_exon_idx_].valid_event(event.type_);
  }

  typename std::vector<size_t>::const_iterator
  begin_introns_for(const Event& event) const {
    const auto& connections_index = introns_for(event);
    return valid_event(event)
      ? connections_index.begin_exon(event.ref_exon_idx_)
      : connections_index.end();
  }
  typename std::vector<size_t>::const_iterator
  end_introns_for(const Event& event) const {
    const auto& connections_index = introns_for(event);
    return valid_event(event)
      ? connections_index.end_exon(event.ref_exon_idx_)
      : connections_index.end();
  }
  typename std::vector<size_t>::const_iterator
  begin_junctions_for(const Event& event) const {
    const auto& connections_index = junctions_for(event);
    return valid_event(event)
      ? connections_index.begin_exon(event.ref_exon_idx_)
      : connections_index.end();
  }
  typename std::vector<size_t>::const_iterator
  end_junctions_for(const Event& event) const {
    const auto& connections_index = junctions_for(event);
    return valid_event(event)
      ? connections_index.end_exon(event.ref_exon_idx_)
      : connections_index.end();
    return junctions_for(event).end_exon(event.ref_exon_idx_);
  }

  bool has_intron(const Event& event) const {
    for (auto it = begin_introns_for(event);
        it != end_introns_for(event); ++it) {
      if ((*introns_)[*it].for_event()) { return true; }
    }
    return false;
  }
  size_t event_size(const Event& event) const {
    size_t ct = 0;
    for (auto it = begin_junctions_for(event);
        it != end_junctions_for(event); ++it) {
      if ((*junctions_)[*it].for_event()) { ++ct; }
    }
    for (auto it = begin_introns_for(event);
        it != end_introns_for(event); ++it) {
      if ((*introns_)[*it].for_event()) { ++ct; }
    }
    return ct;
  }
  bool passed(const Event& event) const {
    for (auto it = begin_junctions_for(event);
        it != end_junctions_for(event); ++it) {
      if ((*junctions_)[*it].for_passed()) { return true; }
    }
    for (auto it = begin_introns_for(event);
        it != end_introns_for(event); ++it) {
      if ((*introns_)[*it].for_passed()) { return true; }
    }
    return false;
  }
  template <bool include_intron>
  std::set<size_t> other_exon_idx_set(
      const Event& event) const {
    std::set<size_t> result;
    for (auto it = begin_junctions_for(event);
        it != end_junctions_for(event); ++it) {
      const auto& x = (*junctions_)[*it];
      if (x.for_event()) { result.insert(x.other_exon_idx(event.type_)); }
    }
    if constexpr(include_intron) {
      for (auto it = begin_introns_for(event);
          it != end_introns_for(event); ++it) {
        const auto& x = (*introns_)[*it];
        if (x.for_event()) { result.insert(x.other_exon_idx(event.type_)); }
      }
    }
    return result;
  }
  /**
   * An event is redundant if it only connects to a single exon
   */
  bool redundant(const Event& event) const {
    constexpr bool INCLUDE_INTRON = true;  // introns count for redundancy
    std::set<size_t> other = other_exon_idx_set<INCLUDE_INTRON>(event);
    return other.size() <= 1;
  }
  /**
   * For strict LSVs, when not redundant or source and mutually redundant
   *
   * This is actually permissive compared to !redundant, because it includes
   * mutually redundant source events
   */
  bool strict_nonredundant(const Event& event) const {
    constexpr bool INCLUDE_INTRON = true;  // introns count for redundancy
    std::set<size_t> other = other_exon_idx_set<INCLUDE_INTRON>(event);
    switch (other.size()) {
      case 0:
        return false;
      case 1:
        // if source event, nonredundant if mutually redundant
        return event.type_ == EventType::SRC_EVENT
          ? redundant(Event{*(other.begin()), EventType::DST_EVENT}) : false;
      default:  // other.size() > 1
        return true;
    }
  }
  /**
   * For permissive LSVs, whenever is source or if target and not redundant
   */
  bool permissive_nonredundant(const Event& event) const {
    if (event.type_ == EventType::SRC_EVENT) {
      return true;
    } else {
      constexpr bool INCLUDE_INTRON = true;  // introns count for redundancy
      std::set<size_t> other = other_exon_idx_set<INCLUDE_INTRON>(event);
      switch (other.size()) {
        case 0:
          return false;
        case 1:
          // this event is redundant, but only false if other event is identical
          return !redundant(Event{*(other.begin()), EventType::SRC_EVENT});
        default:
          return true;
      }
    }
  }
  bool strict_event(const Event& event) const {
    return passed(event) && strict_nonredundant(event);
  }
  bool is_constitutive(const Event& event) const {
    return (event_size(event) == 1) && strict_event(event);
  }
  bool is_strict_LSV(const Event& event) const {
    return (event_size(event) > 1) && strict_event(event);
  }
  bool is_permissive_LSV(const Event& event) const {
    return (event_size(event) > 1) && passed(event)
      && permissive_nonredundant(event);
  }
  bool is_source_LSV(const Event& event) const {
    return event.type_ == EventType::SRC_EVENT
      && (event_size(event) > 1) && passed(event);
  }
  bool is_target_LSV(const Event& event) const {
    return event.type_ == EventType::DST_EVENT
      && (event_size(event) > 1) && passed(event);
  }
  std::string id(const Event& event) const {
    return (*exons_)[event.ref_exon_idx_].event_id(event.type_);
  }
  // non-simplified junction coordinates from event reference exon in direction
  std::vector<position_t> ref_splicesites(const Event& event) const {
    std::vector<position_t> result;
    for (auto it = begin_junctions_for(event);
        it != end_junctions_for(event); ++it) {
      const auto& x = (*junctions_)[*it];
      if (x.simplified()) { continue; }
      const position_t& ss = x.ref_coordinate(event.type_);
      if (result.empty() || result.back() != ss) { result.push_back(ss); }
    }
    return result;
  }
  std::string aborted_description(const Event& event) const {
    std::ostringstream oss;
    oss << event.type_ << (has_intron(event) ? "|na|i" : "|na");
    return oss.str();
  }
  std::string description(const Event& event) const {
    // description =
    // 1. {type()} +
    // 2. |{ref_jidx}e{ref_exct}.{other_pos}o{other_total} for each junction
    // 3. |i if has intron
    // ref_jidx, other_pos ~ enumeration of splicesites on the exons
    // ref_exct ~ how many other exons, other_total ~ splicesites on the other
    // exon
    const std::vector<position_t> ref_ss = ref_splicesites(event);
    std::vector<std::pair<size_t, std::vector<position_t>>> other_exons_ss;
    {  // populate other_exons_ss
      // construct other_exons_ss with events for other exons in other direction
      // TODO(jaicher): verify that INCLUDE_INTRON = false isn't a bug -- I
      // think this would lead to misleading graphics in voila when there are no
      // junctions to the subsequent exon (that the intron is connected to)
      constexpr bool INCLUDE_INTRON = false;  // this is only for junctions
      const std::set<size_t> other_exons
        = other_exon_idx_set<INCLUDE_INTRON>(event);
      const EventType other_type = OtherEventType(event.type_);
      std::transform(other_exons.begin(), other_exons.end(),
          std::back_inserter(other_exons_ss),
          [this, other_type](size_t other_exon_idx) {
          return std::make_pair(
              other_exon_idx,
              ref_splicesites(Event{other_exon_idx, other_type}));
          });
    }  // done constructing other_exons_ss
    // index of splicesites on an exon
    const bool is_forward
      = (*exons_)[event.ref_exon_idx_].gene.strand() == GeneStrandness::FORWARD;
    auto get_ss_idx
      = [is_forward](
          const std::vector<position_t>& ss_vec, position_t query) -> size_t {
        // find splicesite noting in sorted strand order
        auto ss_it = std::lower_bound(ss_vec.begin(), ss_vec.end(), query,
            [is_forward](position_t x, position_t y) {
            return is_forward ? x < y : y < x; });
        return 1 + (ss_it - ss_vec.begin()); };
    // details for another exon
    struct OtherExonDetails { size_t ref_exct, other_pos, other_total; };
    auto get_other_details = [is_forward, &other_exons_ss, &get_ss_idx](
        size_t other_exon_idx, position_t other_position) -> OtherExonDetails {
      OtherExonDetails result;
      auto other_it = std::lower_bound(
          other_exons_ss.begin(), other_exons_ss.end(), other_exon_idx,
          [](const std::pair<size_t, std::vector<position_t>>& x, size_t y) {
          return x.first < y; });
      result.ref_exct = (
          is_forward
          ? 1 + other_it - other_exons_ss.begin()
          : other_exons_ss.end() - other_it);
      const std::vector<position_t>& other_ss_vec = other_it->second;
      result.other_total = other_ss_vec.size();
      result.other_pos = get_ss_idx(other_ss_vec, other_position);
      return result; };
    // create id over junctions/introns that are part of event
    std::ostringstream oss;
    oss << event.type_;
    for (auto it = begin_junctions_for(event);
        it != end_junctions_for(event); ++it) {
      const auto& x = (*junctions_)[*it];
      if (x.for_event()) {
        const auto& [ref_position, other_position]
          = x.coordinates_ref_other(event.type_);
        const size_t& other_exon_idx = x.other_exon_idx(event.type_);
        const size_t ref_jidx = get_ss_idx(ref_ss, ref_position);
        OtherExonDetails other_details
          = get_other_details(other_exon_idx, other_position);
        oss
          << '|' << ref_jidx
          << 'e' << other_details.ref_exct
          << '.' << other_details.other_pos
          << 'o' << other_details.other_total;
        if (oss.tellp() >= EVENT_DESCRIPTION_WIDTH) {
          return aborted_description(event);
        }
      }
    }
    if (has_intron(event)) { oss << "|i"; }
    if (oss.tellp() >= EVENT_DESCRIPTION_WIDTH) {
      return aborted_description(event);
    }
    return oss.str();
  }


  // how do we create Events?
  Events CreateEvents(std::vector<Event>&& events) const {
    std::vector<size_t> connection_offsets{0};
    std::vector<ConnectionIndex> connections;
    for (const Event& event : events) {
      for (auto it = begin_junctions_for(event);
          it != end_junctions_for(event); ++it) {
        if ((*junctions_)[*it].for_event()) {
          constexpr bool IS_INTRON = false;
          connections.push_back(ConnectionIndex{IS_INTRON, *it});
        }
      }
      for (auto it = begin_introns_for(event);
          it != end_introns_for(event); ++it) {
        if ((*introns_)[*it].for_event()) {
          constexpr bool IS_INTRON = true;
          connections.push_back(ConnectionIndex{IS_INTRON, *it});
        }
      }
      connection_offsets.push_back(connections.size());
    }
    return Events{introns_, junctions_, std::move(events),
        std::move(connection_offsets), std::move(connections)};
  }
  Events CreateEvents(std::vector<Event> events) {
    return CreateEvents(std::move(events));
  }
  Events StrictLSVs() const {
    std::vector<Event> events;
    constexpr std::array<EventType, 2> TYPES
      = {EventType::SRC_EVENT, EventType::DST_EVENT};
    for (size_t exon_idx = 0; exon_idx < num_exons(); ++exon_idx) {
      for (const auto& type : TYPES) {
        Event event{exon_idx, type};
        if (is_strict_LSV(event)) { events.push_back(event); }
      }  // loop over event types for a reference exon
    }  // loop over reference exons
    return CreateEvents(std::move(events));
  }
  Events PermissiveLSVs() const {
    std::vector<Event> events;
    constexpr std::array<EventType, 2> TYPES
      = {EventType::SRC_EVENT, EventType::DST_EVENT};
    for (size_t exon_idx = 0; exon_idx < num_exons(); ++exon_idx) {
      for (const auto& type : TYPES) {
        Event event{exon_idx, type};
        if (is_permissive_LSV(event)) { events.push_back(event); }
      }  // loop over event types for a reference exon
    }  // loop over reference exons
    return CreateEvents(std::move(events));
  }
  Events SourceLSVs() const {
    std::vector<Event> events;
    constexpr std::array<EventType, 2> TYPES
      = {EventType::SRC_EVENT, EventType::DST_EVENT};
    for (size_t exon_idx = 0; exon_idx < num_exons(); ++exon_idx) {
      for (const auto& type : TYPES) {
        Event event{exon_idx, type};
        if (is_source_LSV(event)) { events.push_back(event); }
      }  // loop over event types for a reference exon
    }  // loop over reference exons
    return CreateEvents(std::move(events));
  }
  Events TargetLSVs() const {
    std::vector<Event> events;
    constexpr std::array<EventType, 2> TYPES
      = {EventType::SRC_EVENT, EventType::DST_EVENT};
    for (size_t exon_idx = 0; exon_idx < num_exons(); ++exon_idx) {
      for (const auto& type : TYPES) {
        Event event{exon_idx, type};
        if (is_target_LSV(event)) { events.push_back(event); }
      }  // loop over event types for a reference exon
    }  // loop over reference exons
    return CreateEvents(std::move(events));
  }
  Events ConstitutiveEvents() const {
    std::vector<Event> events;
    constexpr std::array<EventType, 2> TYPES
      = {EventType::SRC_EVENT, EventType::DST_EVENT};
    for (size_t exon_idx = 0; exon_idx < num_exons(); ++exon_idx) {
      for (const auto& type : TYPES) {
        Event event{exon_idx, type};
        if (is_constitutive(event)) { events.push_back(event); }
      }  // loop over event types for a reference exon
    }  // loop over reference exons
    return CreateEvents(std::move(events));
  }

 public:
  const std::shared_ptr<Exons>& exons() const { return exons_; }
  const std::shared_ptr<GeneIntrons>& introns() const { return introns_; }
  const std::shared_ptr<GeneJunctions>& junctions() const { return junctions_; }

  const detail::ExonConnectionsIndexes& src_introns() const { return src_introns_; }
  const detail::ExonConnectionsIndexes& dst_introns() const { return dst_introns_; }
  const detail::ExonConnectionsIndexes& src_junctions() const { return src_junctions_; }
  const detail::ExonConnectionsIndexes& dst_junctions() const { return dst_junctions_; }

  ExonConnections(
      const std::shared_ptr<Exons>& exons,
      const std::shared_ptr<GeneIntrons>& introns,
      const std::shared_ptr<GeneJunctions>& junctions)
      : exons_{exons},
        introns_{introns},
        junctions_{junctions},
        src_junctions_{
          detail::ExonConnectionsIndexes::FromConnectionsAndExons(
              *junctions_, exons_, EventType::SRC_EVENT)},
        dst_junctions_{
          detail::ExonConnectionsIndexes::FromConnectionsAndExons(
              *junctions_, exons_, EventType::DST_EVENT)},
        src_introns_{
          detail::ExonConnectionsIndexes::FromConnectionsAndExons(
              *introns_, exons_, EventType::SRC_EVENT)},
        dst_introns_{
          detail::ExonConnectionsIndexes::FromConnectionsAndExons(
              *introns_, exons_, EventType::DST_EVENT)} { }
};

}  // namespace majiq

#endif  // MAJIQ_EXONCONNECTIONS_HPP
