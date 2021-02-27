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
#include <string>
#include <map>
#include <sstream>

#include "Exons.hpp"
#include "GeneJunctions.hpp"
#include "GeneIntrons.hpp"
#include "MajiqConstants.hpp"

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

  std::string id(const Exons& exons) const {
    const Exon& ref_exon = exons[ref_exon_idx_];
    std::ostringstream oss;
    oss << ref_exon.gene.get() << ':'
      << static_cast<char>(event_type_)
      << ':' << ref_exon.coordinates;
    return oss.str();
  }
};
inline bool operator<(const EventIndex& x, const EventIndex& y) noexcept {
  return std::tie(x.ref_exon_idx_, x.event_type_)
    < std::tie(y.ref_exon_idx_, y.event_type_);
}
inline bool operator>(const EventIndex& x, const EventIndex& y) noexcept {
  return y < x;
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

  const std::shared_ptr<Exons>& exons() const { return exons_; }

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
  // if we know for certain that ec_idx is junction or intron
  const GeneJunction& junction(size_t ec_idx) const {
    return (*junctions_)[(*this)[ec_idx].connection_idx_];
  }
  const GeneIntron& intron(size_t ec_idx) const {
    return (*introns_)[(*this)[ec_idx].connection_idx_];
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

  const std::shared_ptr<Exons>& exons() const {
    return event_connections_->exons();
  }

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

  using const_iterator = typename std::vector<EventIndex>::const_iterator;
  const_iterator begin() const { return events_.begin(); }
  const_iterator end() const { return events_.end(); }
  const_iterator find(
      const EventIndex& x, const_iterator lb, const_iterator ub) const {
    // find within known range
    auto it = std::lower_bound(lb, ub, x);
    return (it == end() || *it != x) ? end() : it;
  }
  const_iterator find(const EventIndex& x) const {
    return find(x, begin(), end());
  }
  const_iterator find(
      const EventIndex& x, const_iterator close_it) const {
    if (close_it == end()) { return find(x); }
    const EventIndex& close = *close_it;
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

  const Exon& ref_exon(size_t event_idx) const {
    const size_t exon_idx = events_[event_idx].ref_exon_idx_;
    return (*exons())[exon_idx];
  }
  const KnownGene& gene(size_t event_idx) const {
    return ref_exon(event_idx).gene;
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
  // does the event have any passed connections, regardless of redundancy
  bool passed(size_t event_idx) const {
    for (size_t ec_idx = offsets_[event_idx];
        ec_idx < offsets_[1 + event_idx];
        ++ec_idx) {
      if (event_connections_->for_passed(ec_idx)) { return true; }
    }
    return false;
  }
  std::set<size_t> other_exon_idx_set(
      size_t event_idx, bool include_intron) const {
    std::set<size_t> result;
    for (size_t ec_idx = offsets_[event_idx];
        ec_idx < offsets_[1 + event_idx];
        ++ec_idx) {
      if ((include_intron || !(*event_connections_)[ec_idx].is_intron_)
          && event_connections_->for_event(ec_idx)) {
        result.insert(event_connections_->other_exon_idx(ec_idx));
      }
    }
    return result;
  }
  // is the event redundant?
  bool redundant(size_t event_idx) const {
    constexpr bool INCLUDE_INTRON = true;  // introns matter for redundancy
    std::set<size_t> other = other_exon_idx_set(event_idx, INCLUDE_INTRON);
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
      const auto other_event_it = find(other_event, begin() + event_idx);
      return !(other_event_it == end() || redundant(other_event_it - begin()));
    }
  }
  bool valid_event(size_t event_idx) const {
    return passed(event_idx) && !redundant(event_idx);
  }
  bool is_LSV(size_t event_idx) const {
    return (event_size(event_idx) > 1) && valid_event(event_idx);
  }
  bool is_constitutive(size_t event_idx) const {
    return (event_size(event_idx) == 1) && valid_event(event_idx);
  }
  std::string event_id(size_t event_idx) const {
    return events_[event_idx].id(*exons());
  }
  // NOTE assumes event connections have ref splicesites in strand sorted order
  std::vector<position_t> ref_splicesites(const_iterator event_it) const {
    const size_t event_idx = event_it - begin();
    const EventIndex& event = *event_it;
    // do we use start/end for splice sites?
    const bool strand_forward
      = gene(event_idx).strand() == GeneStrandness::FORWARD;
    const bool is_src = event.event_type_ == EventType::SRC_EVENT;
    const bool use_start = strand_forward == is_src;
    // construct result
    std::vector<position_t> result;
    for (size_t ec_idx = offsets_[event_idx]; ec_idx < offsets_[1 + event_idx];
        ++ec_idx) {
      auto v_connection = event_connections_->junction_or_intron(ec_idx);
      if (std::holds_alternative<GeneJunction>(v_connection)) {
        const auto& j = std::get<GeneJunction>(v_connection);
        if (!j.simplified()) {
          const position_t ss
            = use_start ? j.coordinates.start : j.coordinates.end;
          if (result.empty() || result.back() != ss) {
            result.push_back(ss);
          }  // only if new splice site (using sorted order to only check last)
        }  // only for non-simplified junctions
      }  // only for junctions
    }  // loop over connections
    return result;
  }
  // if we're lazy (or if event_description is too long
  std::string event_description_abort(const size_t event_idx) const {
    std::ostringstream oss;
    oss << static_cast<char>((*this)[event_idx].event_type_)
      << "|na";
    if (has_intron(event_idx)) {
      oss << "|i";
    }
    return oss.str();
  }
  std::string event_description(const size_t event_idx) const {
    // event description = event_type |
    // {ref_jidx}e{ref_exct}.{other_pos}o{other_total} for each junction | i if
    // there is a junction
    // where ref_jidx indexes unique reference splice sites in strand order
    // (including for exitrons, but not simplified)
    // where ref_exct indexes unique other exons (not including self from
    // exitrons) in strand order
    // where other_pos indexes unique splice sites on the other exon (including
    // for exitrons)
    // where other_total tells us how many such splice sites exist on the other
    // exon
    const std::vector<position_t> ref_ss = ref_splicesites(begin() + event_idx);
    std::vector<std::pair<size_t, std::vector<position_t>>> other_exons_ss;
    {
      // construct other_exons_ss with events for other exons in other direction
      constexpr bool INCLUDE_INTRON = false;  // this is only for junctions
      const std::set<size_t> other_exons
        = other_exon_idx_set(event_idx, INCLUDE_INTRON);
      const EventType other_type
        = OtherEventType((*this)[event_idx].event_type_);
      std::transform(other_exons.begin(), other_exons.end(),
          std::back_inserter(other_exons_ss),
          [this, other_type, ref_it = begin() + event_idx](
              size_t other_exon_idx)
              -> std::pair<size_t, std::vector<position_t>> {
          EventIndex other_event = {other_exon_idx, other_type};
          auto other_event_it = find(other_event, ref_it);
          return std::make_pair(
              other_exon_idx, ref_splicesites(other_event_it));
          });
    }  // done constructing other_exons_ss
    // how to get index of a splice site from one of these vectors
    const bool is_forward = gene(event_idx).strand() == GeneStrandness::FORWARD;
    auto get_ss_idx = [is_forward](
        const std::vector<position_t>& ss_vec, position_t query) -> size_t {
      // find splice site noting that it's sorted in strand order
      auto ss_it = std::lower_bound(ss_vec.begin(), ss_vec.end(), query,
          [is_forward](position_t x, position_t y) {
          return is_forward ? x < y : y < x; });
      // index for this purpose starts at 1
      return 1 + (ss_it - ss_vec.begin());
    };
    // how to get details about an other exon
    struct OtherExonDetails {
      size_t ref_exct;
      size_t other_pos;
      size_t other_total;
    };
    auto get_other_details = [is_forward, &other_exons_ss, &get_ss_idx](
        size_t other_exon_idx, position_t other_position) -> OtherExonDetails {
      OtherExonDetails result;
      // get iterator to matching value
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
      return result;
    };
    const bool is_source
      = (*this)[event_idx].event_type_ == EventType::SRC_EVENT;
    // start building description
    std::ostringstream oss;
    oss << static_cast<char>((*this)[event_idx].event_type_);
    for (size_t ec_idx = offsets_[event_idx];
        ec_idx < offsets_[1 + event_idx];
        ++ec_idx) {
      if (event_connections_->for_event(ec_idx)) {
        if ((*event_connections_)[ec_idx].is_intron_) {
          oss << "|i";
        } else {
          const GeneJunction& j = event_connections_->junction(ec_idx);
          position_t ref_position = (
              is_forward == is_source
              ? j.coordinates.start
              : j.coordinates.end);
          position_t other_position = (
              is_forward == is_source
              ? j.coordinates.end
              : j.coordinates.start);
          size_t other_exon_idx = (
              is_source ? j.dst_exon_idx() : j.src_exon_idx());
          const size_t ref_jidx = get_ss_idx(ref_ss, ref_position);
          OtherExonDetails other_details = get_other_details(
              other_exon_idx, other_position);
          oss
            << '|' << ref_jidx
            << 'e' << other_details.ref_exct
            << '.' << other_details.other_pos
            << 'o' << other_details.other_total;
        }
        // abort if too long
        if (oss.tellp() >= EVENT_DESCRIPTION_WIDTH) {
          return event_description_abort(event_idx);
        }
      }  // only consider connections that are for event
    }  // done looping
    return oss.str();
  }

  // TODO(jaicher): simplifier
  // TODO(jaicher): assign coverage to events

  // NOTE(jaicher): do not filter out events. Redundant, non-passed,
  // constitutive events, etc. are necessary for event_description. The one
  // thing we could filter out are simplified event connections
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
