/**
 * EventsAlign.hpp
 *
 * Identify matching events from two splicegraphs
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_EVENTSALIGN_HPP
#define MAJIQ_EVENTSALIGN_HPP

#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Events.hpp"
#include "Exons.hpp"

namespace majiq {

class EventsAlign {
 public:
  struct EventAligned {
    size_t left_idx_;
    size_t right_idx_;
  };
  // events that match from left and right arguments to constructor
  // NOTE: no guarantee of sorting
  std::vector<EventAligned> matched_;

 private:
  /**
   * Check that events match, assuming that known that genes/event type match
   */
  static bool EventsMatch(
      const Events& left, const Events& right,
      size_t left_idx, size_t right_idx) {
    size_t left_start = left.connection_offsets()[left_idx];
    size_t left_size = left.connection_offsets()[1 + left_idx] - left_start;
    size_t right_start = right.connection_offsets()[right_idx];
    {
      size_t right_size
        = right.connection_offsets()[1 + right_idx] - right_start;
      if (right_size != left_size) {
        return false;
      }
    }
    for (size_t i = 0; i < left_size; ++i) {
      if ((left.is_intron(left_start + i) != right.is_intron(right_start + i))
          || (left.connection_start(left_start + i)
            != right.connection_start(right_start + i))
          || (left.connection_end(left_start + i)
            != right.connection_end(right_start + i))) {
        return false;
      }
    }  // done checking for mismatch in each connection
    return true;  // all connections matched, so events match
  }

  template <
    bool IS_INTRON,
    typename ConnectionT
      = std::conditional_t<IS_INTRON, GeneIntron, GeneJunction>,
    typename Compare = detail::CompareContigUnstranded<ConnectionT>>
  static void MatchByConnections(const Events& left, const Events& right,
      std::vector<std::optional<size_t>>& opt_right_idx,
      std::vector<bool>& left_visited,
      std::vector<bool>& right_visited) {
    auto right_idx_it = right.connection_idx_begin<IS_INTRON>();
    const auto right_idx_end = right.connection_idx_end<IS_INTRON>();
    const auto left_idx_end = left.connection_idx_end<IS_INTRON>();
    for (auto left_idx_it = left.connection_idx_begin<IS_INTRON>();
        left_idx_it != left_idx_end; ++left_idx_it) {
      size_t left_eidx = left.connection_event_idx()[*left_idx_it];
      if (left_visited[left_eidx]) {
        continue;  // already has answer
      }
      const auto left_c = left.connection_at<IS_INTRON>(*left_idx_it);
      // get first connection on right that is not less than current on left
      right_idx_it = std::find_if(right_idx_it, right_idx_end,
          [&right, &left_c](const size_t& right_idx) {
          return !Compare{}(
              right.connection_at<IS_INTRON>(right_idx), left_c); });
      for (auto right_idx_eq = right_idx_it;
          (right_idx_eq != right_idx_end
           // while left is also not less than right
           && !Compare{}(
             left_c, right.connection_at<IS_INTRON>(*right_idx_eq)));
          ++right_idx_eq) {
        if (left_c.gene == right.connection_at<IS_INTRON>(*right_idx_eq).gene
            && (left.connection_event(*left_idx_it).type_
              == right.connection_event(*right_idx_eq).type_)) {
          // same event-connection! check...
          left_visited[left_eidx] = true;
          size_t right_eidx = right.connection_event_idx()[*right_idx_it];
          if (!right_visited[right_eidx]) {
            right_visited[right_eidx] = true;
            if (EventsMatch(left, right, left_eidx, right_eidx)) {
              opt_right_idx[left_eidx] = right_eidx;
            }  // update index from null if events match
          }  // only can match if right event hasn't been visited before
          break;  // found unique potential match for event connection
        }
      }  // loop over right indexes that maatch contig and coordinates
    }  // loop over left indexes in contig region sorted order
  }

 public:
  EventsAlign(const Events& left, const Events& right) : matched_{} {
    if (left.num_events() > right.num_events()) {
      auto reversed = EventsAlign{right, left};
      matched_ = std::move(reversed.matched_);
      for (auto&& x : matched_) {
        std::swap(x.left_idx_, x.right_idx_);
      }
      return;
    }
    // otherwise
    if (left.introns()->parents() != right.introns()->parents()) {
      throw std::runtime_error("left and right events must share genes");
    }
    // indexed along left events, indicate matching index in right events
    std::vector<std::optional<size_t>> opt_right_idx(left.num_events());
    std::vector<bool> left_visited(left.num_events());
    std::vector<bool> right_visited(right.num_events());
    // update these using introns and junctions
    MatchByConnections<true>(
        left, right, opt_right_idx, left_visited, right_visited);
    MatchByConnections<false>(
        left, right, opt_right_idx, left_visited, right_visited);
    // create matched_
    for (size_t left_eidx = 0; left_eidx < left.num_events(); ++left_eidx) {
      if (opt_right_idx[left_eidx].has_value()) {
        matched_.push_back(
            EventAligned{left_eidx, *(opt_right_idx[left_eidx])});
      }
    }
    return;
  }
};  // struct EventsAlign

}  // namespace majiq


#endif  // MAJIQ_EVENTSALIGN_HPP
