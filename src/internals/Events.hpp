/**
 * Events.hpp
 *
 * Events that reflect a source or target exon shared by junctions and/or
 * intron
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_EVENTS_HPP
#define MAJIQ_EVENTS_HPP

#include <algorithm>
#include <memory>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "MajiqTypes.hpp"
#include "GeneIntrons.hpp"
#include "GeneJunctions.hpp"

namespace majiq {

struct ConnectionIndex {
  bool is_intron_;
  size_t idx_;
};

class Events {
 private:
  // underlying connections
  const std::shared_ptr<GeneIntrons> introns_;
  const std::shared_ptr<GeneJunctions> junctions_;

  // indexes for events, for connections that are part of each event
  const std::vector<Event> events_;
  const std::vector<size_t> connection_offsets_;
  const std::vector<ConnectionIndex> connections_;

  // get indexes of connections_ into junctions/introns, sorted by contig/region
  const std::vector<size_t> junction_connection_idx_;
  const std::vector<size_t> intron_connection_idx_;

 public:
  size_t num_events() const { return events_.size(); }
  size_t num_connections() const { return connections_.size(); }
  size_t num_junctions() const { return junction_connection_idx_.size(); }
  size_t num_introns() const { return intron_connection_idx_.size(); }
  size_t size() const { return num_connections(); }

  // index into connections for contig-sorted junctins
  const size_t& sorted_junction_connection_idx(size_t sorted_idx) const {
    return junction_connection_idx_[sorted_idx];
  }
  const size_t& sorted_intron_connection_idx(size_t sorted_idx) const {
    return intron_connection_idx_[sorted_idx];
  }

  // get underlying junction or intron from connection_idx
  // NOTE: does not check that actually a junction or intron
  const GeneJunction& connection_junction(size_t connection_idx) const {
    return (*junctions_)[connections_[connection_idx].idx_];
  }
  const GeneIntron& connection_intron(size_t connection_idx) const {
    return (*introns_)[connections_[connection_idx].idx_];
  }
  // check if intron or not
  const bool& is_intron(size_t connection_idx) const {
    return connections_[connection_idx].is_intron_;
  }

  const std::vector<Event>& events() const { return events_; }
  const std::vector<size_t>& connection_offsets() const {
    return connection_offsets_;
  }
  const std::vector<ConnectionIndex>& connections() const {
    return connections_;
  }

 private:
  // either for introns or junctions
  template <bool INTRON>
  std::vector<size_t> ContigSortedConnectionIndexes() {
    using ConnectionT = std::conditional_t<INTRON, GeneIntron, GeneJunction>;
    auto get_connection = [this](size_t connection_idx) -> const ConnectionT& {
      if constexpr(INTRON) {
        return connection_intron(connection_idx);
      } else {
        return connection_junction(connection_idx);
      } };
    // get connection_idx matching is_intron
    std::vector<size_t> result;
    for (size_t i = 0; i < num_connections(); ++i) {
      if (is_intron(i) == INTRON) { result.push_back(i); }
    }
    // sort result as unstranded contig region
    std::sort(result.begin(), result.end(),
        [&get_connection](size_t i, size_t j) -> bool {
        return detail::CompareContigUnstranded<ConnectionT>{}(
            get_connection(i), get_connection(j)); });
    // return result
    return result;
  }

 public:
  Events(const std::shared_ptr<GeneIntrons>& introns,
      const std::shared_ptr<GeneJunctions>& junctions,
      std::vector<Event>&& events,
      std::vector<size_t>&& connection_offsets,
      std::vector<ConnectionIndex>&& connections)
      : introns_{introns},
        junctions_{junctions},
        events_{std::move(events)},
        connection_offsets_{std::move(connection_offsets)},
        connections_{std::move(connections)},
        junction_connection_idx_{ContigSortedConnectionIndexes<false>()},
        intron_connection_idx_{ContigSortedConnectionIndexes<true>()} { }
};


};

#endif  // MAJIQ_EVENTS_HPP
