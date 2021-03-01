/**
 * Events.hpp
 *
 * Container that manages events and event connections together
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_EVENTS_HPP
#define MAJIQ_EVENTS_HPP

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "EventConnection.hpp"
#include "EventConnections.hpp"


namespace majiq {

class EventReference;  // defined after Events

class Events {
 private:
  const std::shared_ptr<EventConnections> event_connections_;
  const std::vector<size_t> offsets_;  // into event connections
  const std::vector<Event> events_;
  const std::vector<size_t> parent_idx_offsets_;  // genes into this

 public:
  const std::shared_ptr<EventConnections>& event_connections() const {
    return event_connections_;
  }
  const std::vector<size_t>& offsets() const { return offsets_; }
  const std::vector<Event>& events() const { return events_; }
  const std::shared_ptr<Genes>& genes() const {
    return event_connections_->exons()->parents();
  }
  const std::shared_ptr<Genes>& parents() const { return genes(); }
  const std::vector<size_t>& parent_idx_offsets() const {
    return parent_idx_offsets_;
  }
  const std::vector<size_t>& gene_idx_offsets() const {
    return parent_idx_offsets();
  }

  size_t size() const { return events_.size(); }

  using const_iterator = typename std::vector<Event>::const_iterator;
  EventReference operator[](size_t idx) const;
  EventReference operator[](const_iterator event_it) const;

  const_iterator begin() const { return events_.cbegin(); }
  const_iterator end() const { return events_.cend(); }
  const_iterator begin_parent(size_t gene_idx) const {
    return begin() + parent_idx_offsets_[gene_idx];
  }
  const_iterator end_parent(size_t gene_idx) const {
    return begin_parent(1 + gene_idx);
  }
  const_iterator begin_parent(const KnownGene& gene) const {
    return begin_parent(gene.idx_);
  }
  const_iterator end_parent(const KnownGene& gene) const {
    return end_parent(gene.idx_);
  }

  const_iterator find(
      const Event& x, const_iterator lb, const_iterator ub) const {
    // find within known range
    auto it = std::lower_bound(lb, ub, x);
    return (it == end() || *it != x) ? end() : it;
  }
  const_iterator find(const Event& x) const {
    return find(x, begin(), end());
  }
  const_iterator find(
      const Event& x, const_iterator close_it) const {
    if (close_it == end()) { return find(x); }
    const Event& close = *close_it;
    if (close == x) {
      return close_it;
    }
    // otherwise, define lb, ub to do search
    const_iterator lb, ub;
    if (close < x) {
      lb = 1 + close_it;
      std::ptrdiff_t max_distance = std::min(
          end() - lb,
          static_cast<std::ptrdiff_t>(
            1 + 2 * (x.ref_exon_idx_ - close.ref_exon_idx_)));
      ub = lb + max_distance;
    } else {
      ub = close_it;
      std::ptrdiff_t max_distance = std::min(
          ub - begin(),
          static_cast<std::ptrdiff_t>(
            1 + 2 * (close.ref_exon_idx_ - x.ref_exon_idx_)));
      lb = ub - max_distance;
    }
    return find(x, lb, ub);
  }

  friend class EventReference;

 private:
  static std::vector<size_t> CalculateOffsets(
      const std::vector<EventConnection>& event_connections) {
    std::vector<size_t> result = {0};
    auto prev_event_it = event_connections.begin();
    for (auto event_it = event_connections.begin();
        event_it != event_connections.end();
        ++event_it) {
      if (prev_event_it->event_ != event_it->event_) {
        result.push_back(event_it - event_connections.begin());
        prev_event_it = event_it;
      }
    }
    result.push_back(event_connections.size());
    return result;
  }
  static std::vector<Event> GetIndex(
      const std::vector<EventConnection>& event_connections,
      const std::vector<size_t>& event_offsets) {
    std::vector<Event> result(event_offsets.size() - 1);
    std::transform(event_offsets.begin(), event_offsets.end() - 1,
        result.begin(),
        [&event_connections](size_t i) {
        return event_connections[i].event_; });
    return result;
  }
  static std::vector<size_t> GetGeneOffsets(
      const Exons& exons, const std::vector<Event>& events) {
    std::vector<size_t> result(1 + exons.parents()->size(), 0);
    size_t gene_idx = 0;
    for (auto ev_it = events.begin(); ev_it != events.end(); ++ev_it) {
      const size_t ev_gene_idx = exons[ev_it->ref_exon_idx_].gene.idx_;
      if (ev_gene_idx < gene_idx) {
        throw std::logic_error("Events are not sorted with respect to genes");
      }
      for (; gene_idx < ev_gene_idx; ++gene_idx) {
        result[1 + gene_idx] = ev_it - events.begin();
      }
    }
    for (; gene_idx < exons.parents()->size(); ++gene_idx) {
      result[1 + gene_idx] = events.size();
    }
    return result;
  }

 public:
  Events(const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<GeneIntrons>& introns)
      : event_connections_{
          std::make_shared<EventConnections>(junctions, introns)},
        offsets_{CalculateOffsets(event_connections_->event_connections())},
        events_{GetIndex(event_connections_->event_connections(), offsets_)},
        parent_idx_offsets_{
          GetGeneOffsets(*event_connections_->exons(), events_)} { }

  // TODO(jaicher): simplifier
  // TODO(jaicher): assign coverage to events

  // NOTE(jaicher): do not filter out events. Redundant, non-passed,
  // constitutive events, etc. are necessary for event_description. The one
  // thing we could filter out are simplified event connections
};

class EventReference {
 private:
  const Events& parent_;
  size_t idx_;

  const EventConnections& connections() const {
    return *parent_.event_connections();
  }
  size_t begin_idx() const { return parent_.offsets_[idx_]; }
  size_t end_idx() const { return parent_.offsets_[1 + idx_]; }

 public:
  EventReference(const Events& parent, size_t idx)
      : parent_{parent}, idx_{idx} { }

  size_t idx() const { return idx_; }

  const Event& event() const { return parent_.events_[idx_]; }
  size_t ref_exon_idx() const { return event().ref_exon_idx_; }
  EventType type() const { return event().event_type_; }
  bool is_source() const { return type() == EventType::SRC_EVENT; }
  const Exon& ref_exon() const {
    return (*connections().exons())[ref_exon_idx()];
  }
  const KnownGene& gene() const { return ref_exon().gene; }
  GeneStrandness strand() const { return gene().strand(); }
  bool is_forward() const { return strand() == GeneStrandness::FORWARD; }
  std::string id() const {
    const Exon& x = ref_exon();
    std::ostringstream oss;
    oss << x.gene.get() << ':'
      << static_cast<char>(type())
      << ':' << x.coordinates;
    return oss.str();
  }

  bool has_intron() const { return connections()[end_idx() - 1].is_intron(); }
  size_t event_size() const {
    size_t ct = 0;
    for (size_t i = begin_idx(); i < end_idx(); ++i) {
      if (connections()[i].for_event()) { ++ct; }
    }
    return ct;
  }
  bool passed() const {
    for (size_t i = begin_idx(); i < end_idx(); ++i) {
      if (connections()[i].for_passed()) { return true; }
    }
    return false;
  }
  bool redundant() const;
  bool valid_event() const;
  bool is_LSV() const;
  bool is_constitutive() const;

  // get unique other exon
  std::set<size_t> other_exon_idx_set(bool include_intron) const {
    std::set<size_t> result;
    for (size_t i = begin_idx(); i < end_idx(); ++i) {
      auto x = connections()[i];
      if ((include_intron || x.is_junction()) && x.for_event()) {
        result.insert(x.other_exon_idx());
      }
    }
    return result;
  }
  // get unique reference positions from junctions
  std::vector<position_t> ref_splicesites() const {
    // NOTE: reference positions are in strand order for an event, so we can
    // easily filter duplicates
    const bool use_start = is_source() == is_forward();
    std::vector<position_t> result;
    for (size_t i = begin_idx(); i < end_idx(); ++i) {
      auto x = connections()[i];
      if (x.is_intron()) { continue; }
      const auto& j = x.junction();
      if (j.simplified()) { continue; }
      const position_t ss
        = use_start ? j.coordinates.start : j.coordinates.end;
      if (result.empty() || result.back() != ss) {
        result.push_back(ss);
      }  // if new splicesite
    }  // loop over event connections
    return result;
  }

  std::string aborted_description() const {
    std::ostringstream oss;
    oss << static_cast<char>(type()) << "|na";
    if (has_intron()) { oss << "|i"; }
    return oss.str();
  }
  std::string description() const;
};

inline EventReference Events::operator[](size_t idx) const {
  return EventReference{*this, idx};
}
inline EventReference Events::operator[](Events::const_iterator it) const {
  return (*this)[it - begin()];
}

inline bool EventReference::redundant() const {
  constexpr bool INCLUDE_INTRON = true;  // introns matter for redundancy
  std::set<size_t> other = other_exon_idx_set(INCLUDE_INTRON);
  if (other.size() == 0) {
    // shouldn't happen, but maybe if someone simplified everything?
    return true;
  } else if (other.size() > 1) {
    return false;
  } else if (!is_source()) {
    return true;
  } else {
    // source events with 1 exon can be nonredundant if matching target only
    // has source as its exon
    const auto other_event = Event{*other.begin(), EventType::DST_EVENT};
    const auto other_it = parent_.find(other_event, parent_.begin() + idx());
    return !(other_it == parent_.end() || parent_[other_it].redundant());
  }
}
inline bool EventReference::valid_event() const {
  return passed() && !redundant();
}
inline bool EventReference::is_LSV() const {
  return (event_size() > 1) &&  valid_event();
}
inline bool EventReference::is_constitutive() const {
  return (event_size() == 1) &&  valid_event();
}
inline std::string EventReference::description() const {
  // description =
  // 1. {type()} +
  // 2. |{ref_jidx}e{ref_exct}.{other_pos}o{other_total} for each junction
  // 3. |i if has intron
  // ref_jidx, other_pos ~ enumeration of splicesites on the exons
  // ref_exct ~ how many other exons, other_total ~ splicesites on the other
  // exon
  const std::vector<position_t> ref_ss = ref_splicesites();
  std::vector<std::pair<size_t, std::vector<position_t>>> other_exons_ss;
  {  // populate other_exons_ss
    // construct other_exons_ss with events for other exons in other direction
    // TODO(jaicher): verify that INCLUDE_INTRON = false isn't a bug -- I
    // think this would lead to misleading graphics in voila when there are no
    // junctions to the subsequent exon (that the intron is connected to)
    constexpr bool INCLUDE_INTRON = false;  // this is only for junctions
    const std::set<size_t> other_exons = other_exon_idx_set(INCLUDE_INTRON);
    const EventType other_type = OtherEventType(type());
    std::transform(other_exons.begin(), other_exons.end(),
        std::back_inserter(other_exons_ss),
        [this, other_type, ref_it = parent_.begin() + idx_](
            size_t other_exon_idx)
            -> std::pair<size_t, std::vector<position_t>> {
          Event other_event = {other_exon_idx, other_type};
          auto other_event_it = parent_.find(other_event, ref_it);
          return std::make_pair(
              other_exon_idx, parent_[other_event_it].ref_splicesites());
        });
  }  // done constructing other_exons_ss
  // index of splicesites on an exon
  const bool is_forward = this->is_forward();
  const bool is_source = this->is_source();
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
  auto get_other_details
    = [is_forward, &other_exons_ss, &get_ss_idx](
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
  std::ostringstream oss;
  oss << static_cast<char>(type());
  for (size_t i = begin_idx(); i < end_idx(); ++i) {
    auto x = connections()[i];
    if (!x.for_event()) {
      continue;  // ignore connections not for event (i.e. simplified, exitron)
    }
    if (x.is_intron()) {
      oss << "|i";
    } else {
      const GeneJunction& j = x.junction();
      position_t ref_position, other_position;
      size_t other_exon_idx;
      if (is_forward == is_source) {
        ref_position = j.coordinates.start;
        other_position = j.coordinates.end;
        other_exon_idx = j.end_exon_idx();
      } else {
        ref_position = j.coordinates.end;
        other_position = j.coordinates.start;
        other_exon_idx = j.start_exon_idx();
      }
      const size_t ref_jidx = get_ss_idx(ref_ss, ref_position);
      OtherExonDetails other_details
        = get_other_details(other_exon_idx, other_position);
      oss
        << '|' << ref_jidx
        << 'e' << other_details.ref_exct
        << '.' << other_details.other_pos
        << 'o' << other_details.other_total;
    }  // junction vs intron
    // abort if too long
    if (oss.tellp() >= EVENT_DESCRIPTION_WIDTH) {
      return aborted_description();
    }
  }  // done looping over connections for the event
  return oss.str();
}

}  // namespace majiq

#endif  // MAJIQ_EVENTS_HPP
