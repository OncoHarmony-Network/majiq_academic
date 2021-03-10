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
  typename std::vector<size_t>::const_iterator
  connection_idx_junctions_begin() const {
    return junction_connection_idx_.cbegin();
  }
  typename std::vector<size_t>::const_iterator
  connection_idx_junctions_end() const {
    return junction_connection_idx_.cend();
  }
  typename std::vector<size_t>::const_iterator
  connection_idx_introns_begin() const {
    return intron_connection_idx_.cbegin();
  }
  typename std::vector<size_t>::const_iterator
  connection_idx_introns_end() const {
    return intron_connection_idx_.cend();
  }

  // get underlying junction or intron from connection_idx
  // NOTE: does not check that actually a junction or intron
  const GeneJunction& connection_junction(size_t connection_idx) const {
    return (*junctions_)[connections_[connection_idx].idx_];
  }
  const GeneIntron& connection_intron(size_t connection_idx) const {
    return (*introns_)[connections_[connection_idx].idx_];
  }

  // provide templated access into junctions or introns
  template <bool IS_INTRON>
  const std::conditional_t<IS_INTRON, GeneIntron, GeneJunction>&
  connection_at(size_t connection_idx) const {
    if constexpr(IS_INTRON) {
      return connection_intron(connection_idx);
    } else {
      return connection_junction(connection_idx);
    }
  }
  template <bool IS_INTRON>
  typename std::vector<size_t>::const_iterator
  connection_idx_begin() const {
    if constexpr(IS_INTRON) {
      return connection_idx_introns_begin();
    } else {
      return connection_idx_junctions_begin();
    }
  }
  template <bool IS_INTRON>
  typename std::vector<size_t>::const_iterator
  connection_idx_end() const {
    if constexpr(IS_INTRON) {
      return connection_idx_introns_end();
    } else {
      return connection_idx_junctions_end();
    }
  }

  // check if intron or not
  const bool& is_intron(size_t connection_idx) const {
    return connections_[connection_idx].is_intron_;
  }

  // get region information from connection_idx
  const KnownGene& connection_gene(size_t connection_idx) const {
    return is_intron(connection_idx)
      ? connection_intron(connection_idx).gene
      : connection_junction(connection_idx).gene;
  }
  const position_t& connection_start(size_t connection_idx) const {
    return is_intron(connection_idx)
      ? connection_intron(connection_idx).coordinates.start
      : connection_junction(connection_idx).coordinates.start;
  }
  const position_t& connection_end(size_t connection_idx) const {
    return is_intron(connection_idx)
      ? connection_intron(connection_idx).coordinates.end
      : connection_junction(connection_idx).coordinates.end;
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
        intron_connection_idx_{ContigSortedConnectionIndexes<true>()} {
    if (introns == nullptr) {
      throw std::runtime_error("Events given null introns");
    } else if (junctions == nullptr) {
      throw std::runtime_error("Events given null introns");
    } else if (introns->parents() != junctions->parents()) {
      throw std::runtime_error("Event junctions/introns do not share genes");
    }
  }
};


};

#endif  // MAJIQ_EVENTS_HPP
