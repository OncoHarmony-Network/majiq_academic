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
 public:
  using event_iterator = typename std::vector<Event>::const_iterator;

 private:
  // underlying connections
  const std::shared_ptr<GeneIntrons> introns_;
  const std::shared_ptr<GeneJunctions> junctions_;
  const std::shared_ptr<Exons> exons_;

  // event definitions: ref_exon and type (source vs target)
  const std::vector<Event> events_;

  // go from genes to events
  // require that events_ are sorted and unique, start at 0, end at events_.size
  const std::vector<size_t> gene_idx_offsets_;

  // connection definitions: is intron/junction, index into introns/junctions
  const std::vector<ConnectionIndex> connections_;

  // go from events to connections using offsets
  // require sorted, start/end at 0/connections_.size(), size=1+events_.size
  const std::vector<size_t> connection_offsets_;

  // go from connections back to events they belong to
  const std::vector<size_t> connection_event_idx_;

  // get indexes of connections_ into junctions/introns, sorted by contig/region
  // Used by EventsCoverage to assign coverage in linear sweep over SJ, which
  // are also sorted by contig/region
  const std::vector<size_t> junction_connection_idx_;
  const std::vector<size_t> intron_connection_idx_;

 public:
  event_iterator events_begin() const { return events_.begin(); }
  event_iterator events_end() const { return events_.end(); }
  event_iterator events_begin_gene(size_t gene_idx) const {
    return events_begin() + gene_idx_offsets_[gene_idx];
  }
  event_iterator events_end_gene(size_t gene_idx) const {
    return events_begin_gene(1 + gene_idx);
  }
  /**
   * first iterator over events_ that is not less than query
   */
  event_iterator events_lower_bound(const Event& query) const {
    // verify that valid reference exon
    if (query.ref_exon_idx_ >= exons_->size()) {
      throw std::invalid_argument(
          "events_lower_bound query has invalid ref_exon_idx");
    }
    // get gene associated with event
    const auto gene_idx = (*exons_)[query.ref_exon_idx_].gene.idx_;
    // perform query on the range associated with matching gene
    return std::lower_bound(
        events_begin_gene(gene_idx), events_end_gene(gene_idx), query);
  }


  size_t num_events() const { return events_.size(); }
  size_t num_connections() const { return connections_.size(); }
  size_t num_junctions() const { return junction_connection_idx_.size(); }
  size_t num_introns() const { return intron_connection_idx_.size(); }
  size_t size() const { return num_connections(); }

  const std::shared_ptr<GeneIntrons>& introns() const { return introns_; }
  const std::shared_ptr<GeneJunctions>& junctions() const { return junctions_; }
  const std::shared_ptr<Exons>& exons() const { return exons_; }

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
  const Exon& event_ref_exon(size_t event_idx) const {
    return (*exons_)[events_[event_idx].ref_exon_idx_];
  }
  const Exon& connection_ref_exon(size_t connection_idx) const {
    return event_ref_exon(connection_event_idx_[connection_idx]);
  }
  std::string event_id(size_t event_idx) const {
    return event_ref_exon(event_idx).event_id(events_[event_idx].type_);
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
    return connection_ref_exon(connection_idx).gene;
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
  const bool& connection_denovo(size_t connection_idx) const {
    return is_intron(connection_idx)
      ? connection_intron(connection_idx).denovo()
      : connection_junction(connection_idx).denovo();
  }

  const Event& connection_event(size_t connection_idx) const {
    return events_[connection_event_idx_[connection_idx]];
  }
  const size_t& connection_other_exon_idx(size_t connection_idx) const {
    const EventType& type = connection_event(connection_idx).type_;
    return is_intron(connection_idx)
      ? connection_intron(connection_idx).other_exon_idx(type)
      : connection_junction(connection_idx).other_exon_idx(type);
  }

  const std::vector<Event>& events() const { return events_; }
  const std::vector<size_t>& gene_idx_offsets() const {
    return gene_idx_offsets_;
  }
  const std::vector<size_t>& connection_offsets() const {
    return connection_offsets_;
  }
  const std::vector<size_t>& connection_event_idx() const {
    return connection_event_idx_;
  }
  const std::vector<ConnectionIndex>& connections() const {
    return connections_;
  }

 private:
  /**
   * Get connection indexes corresponding to introns or junctions in unstranded
   * contig sorted order
   *
   * Checks that
   * - introns/junctions are not null
   * - connection indexes are valid into introns/junctions
   * - connection gene matches that of exon
   *
   * Assumes that gene_idx_offsets_ already set with GeneOffsetsForEvents,
   * which checks validity of exon indexes
   * Assumes that connection_event_idx already set with
   * ConnectionEventIndexFromOffsets
   */
  template <bool INTRON>
  std::vector<size_t> ContigSortedConnectionIndexes() const {
    // verify that introns/junctions not null
    if constexpr(INTRON) {
      if (introns_ == nullptr) {
        throw std::invalid_argument("Events given null introns");
      }
    } else {
      if (junctions_ == nullptr) {
        throw std::invalid_argument("Events given null junctions");
      }
    }
    // function to get intron/junction for specified index
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
      // check that the connection references a valid intron/junction
      if (is_intron(i) == INTRON) {
        if constexpr(INTRON) {
          if (connections_[i].idx_ >= introns_->size()) {
            throw std::invalid_argument(
                "Events has invalid index into introns");
          }
        } else {
          if (connections_[i].idx_ >= junctions_->size()) {
            throw std::invalid_argument(
                "Events has invalid index into junctions");
          }
        }
        // check that connection has same gene as corresponding event
        if (get_connection(i).gene != connection_ref_exon(i).gene) {
          throw std::invalid_argument(
              "Events has event/connection with non-matching genes");
        }
        // this index has the right type, so track it
        result.push_back(i);
      }
    }
    // sort result as unstranded contig region
    std::sort(result.begin(), result.end(),
        [&get_connection](size_t i, size_t j) -> bool {
        return detail::CompareContigUnstranded<ConnectionT>{}(
            get_connection(i), get_connection(j)); });
    // return result
    return result;
  }

  /**
   * Map connections back to events using connection offsets (from events)
   *
   * Checks that connection_offsets_ is valid: starts at zero, ends at number
   * of connections, nondecreasing, and has correct length (1 + events.size)
   *
   * Returns vector with same length as connections. Each value corresponds to
   * the event that connection belongs to
   */
  std::vector<size_t> ConnectionEventIndexFromOffsets() const {
    if (connection_offsets_.size() != 1 + num_events()) {
      throw std::invalid_argument(
          "Events connection offsets have incorrect size for events");
    } else if (connection_offsets_[0] != 0) {
      throw std::invalid_argument(
          "Events connection offsets must start at 0");
    } else if (connection_offsets_.back() != num_connections()) {
      throw std::invalid_argument(
          "Events connection offsets must end at number of connections");
    }
    std::vector<size_t> result(num_connections());
    for (size_t i = 0; i < num_events(); ++i) {
      if (connection_offsets_[i] > connection_offsets_[1 + i]
          || connection_offsets_[1 + i] > result.size()) {
        throw std::invalid_argument(
            "Event connection offsets must be nondecreasing");
      }
      std::fill(result.begin() + connection_offsets_[i],
          result.begin() + connection_offsets_[1 + i], i);
    }
    return result;
  }

  /**
   * Check that events are in sorted order, use exons to get offsets for genes
   */
  static std::vector<size_t> GeneOffsetsForEvents(
      const std::vector<Event>& events,
      const std::shared_ptr<Exons>& exons_ptr) {
    // get exons events are defined over
    if (exons_ptr == nullptr) {
      throw std::invalid_argument("Events has null exons");
    }
    const Exons& exons = *exons_ptr;
    const size_t num_genes = exons.parents()->size();

    // offsets over genes must have length 1 + number of genes
    std::vector<size_t> result(1 + num_genes);

    // check that events are sorted and valid, update offsets
    size_t cur_gene_idx = 0;
    for (size_t i = 0; i < events.size(); ++i) {
      const size_t ref_exon_idx = events[i].ref_exon_idx_;
      if (ref_exon_idx >= exons.size()) {
        throw std::invalid_argument("Events given invalid ref_exon_idx");
      }
      const size_t new_gene_idx = exons[ref_exon_idx].gene.idx_;
      // by construction of Exons, cur_gene_idx is valid
      for (; cur_gene_idx < new_gene_idx; ++cur_gene_idx) {
        result[1 + cur_gene_idx] = i;
      }
      // validate that events are sorted
      if (i > 0 && !(events[i - 1] < events[i])) {
        throw std::invalid_argument("Events are not in sorted order");
      }
    }
    // fill in remaining offsets
    for (; cur_gene_idx < num_genes; ++cur_gene_idx) {
      result[1 + cur_gene_idx] = events.size();
    }
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
        exons_{introns_->connected_exons()},
        // event definitions without checking anything
        events_{std::move(events)},
        // check exons not null, events sorted, valid exon indexes
        gene_idx_offsets_{GeneOffsetsForEvents(events_, exons_)},
        // connection/offset definitions without checking anything
        connections_{std::move(connections)},
        connection_offsets_{std::move(connection_offsets)},
        // checks that connection_offsets_ defined appropriately
        connection_event_idx_{ConnectionEventIndexFromOffsets()},
        // introns/junctions not null, connections valid, have matching genes
        junction_connection_idx_{ContigSortedConnectionIndexes<false>()},
        intron_connection_idx_{ContigSortedConnectionIndexes<true>()} {
    // recall that exons_ is defined using introns_->connected_exons()
    // this also checks that they share genes because they can only have
    // connected exons set if they share genes with exons
    if (exons_ != junctions_->connected_exons()) {
      throw std::runtime_error("Events junctions/introns do not share exons");
    }
  }
  Events(const Events&) = default;
  Events(Events&&) = default;
  Events& operator=(const Events&) = default;
  Events& operator=(Events&&) = default;
};


};

#endif  // MAJIQ_EVENTS_HPP
