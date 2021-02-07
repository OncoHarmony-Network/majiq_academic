/**
 * Regions.hpp
 *
 * Regions base class for splicegraph (closed intervals, relative to known
 * contigs or known genes)
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_REGIONS_HPP
#define MAJIQ_REGIONS_HPP

#include <iostream>
#include <tuple>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <memory>
#include <functional>
#include <initializer_list>
#include <utility>

#include "ContigRegion.hpp"
#include "GeneRegion.hpp"


namespace majiq {
namespace detail {

template <typename RegionT, bool HAS_OVERLAPS>
class Regions {
 public:
  using ParentT = decltype(std::declval<RegionT>().parent());
  using ParentsPtrT = decltype(std::declval<ParentT>().ptr_);
  using IntervalT = decltype(std::declval<RegionT>().coordinates);
  using DataT = decltype(std::declval<RegionT>().data);
  using value_type = RegionT;
  using vecT = std::vector<value_type>;
  using const_iterator = typename vecT::const_iterator;

 public:
  const vecT elements_;
  const ParentsPtrT parents_;
  // since sorted, we can preprocess to know where each of parents are
  const std::vector<size_t> parent_idx_offsets_;
  // if there are overlaps we can know the max position seen so ends are
  // ordered too (NOTE: assumes sorted within parent by (start, end))
  const std::vector<position_t> elements_end_cummax_;

 private:
  static std::vector<size_t> validate_and_parent_offsets(
      const vecT& elements, const ParentsPtrT& parents) {
    // parent offsets are 0 if no parents
    if (parents == nullptr) { return {0}; }
    // otherwise, we at least know the length of the result
    std::vector<size_t> result(parents->size());
    // loop through elements, noting when we have a new parent
    // also, validate that it is sorted
    size_t cur_parent_idx = 0;
    for (size_t i = 0; i < elements.size(); ++i) {
      // update result tracking offsets
      const size_t new_parent_idx = elements[i].parent().idx_;
      for (; cur_parent_idx < new_parent_idx; ++cur_parent_idx) {
        result[1 + cur_parent_idx] = i;
      }
      // validate sorted, share same pointer to parents
      if (parents != elements[i].parent().ptr_) {
        throw std::invalid_argument(
            "Regions must all point to same container of parents");
      }
      if (i > 0 && !(elements[i - 1] < elements[i])) {
        throw std::invalid_argument("Regions must be in sorted order");
      }
    }
    // fill in remaining offsets to end of elements
    for (; cur_parent_idx < parents->size(); ++cur_parent_idx) {
      result[1 + cur_parent_idx] = elements.size();
    }
    return result;
  }
  static std::vector<position_t> get_cummax_ends(const vecT& elements) {
    if constexpr(!HAS_OVERLAPS) {
      // no need for this, so just leave it empty
      return std::vector<position_t>{};
    } else {
      // we know the size of our result
      std::vector<position_t> result(elements.size());
      if (result.empty()) { return result; }
      // first element is already known
      result[0] = elements[0].coordinates.last_pos();
      for (size_t i = 1; i < elements.size(); ++i) {
        result[i] = (
            (elements[i - 1].parent() != elements[i].parent())
            || (result[i - 1] < elements[i].coordinates.last_pos()))
          ? elements[i].coordinates.last_pos() : result[i - 1];
      }
      return result;
    }
  }

 public:
  const size_t size() const noexcept { return elements_.size(); }
  bool empty() const noexcept { return elements_.empty(); }
  const value_type& operator[](size_t idx) const { return elements_[idx]; }
  const_iterator begin() const { return elements_.cbegin(); }
  const_iterator end() const { return elements_.cend(); }
  const ParentsPtrT& parents() const { return parents_; }


  const_iterator find(const RegionT& key) const {
    // iterator into elements_ that are first vs last
    const_iterator first = begin() + parent_idx_offsets_[key.parent().idx_];
    const_iterator last = begin() + parent_idx_offsets_[1 + key.parent().idx_];
    // do search on this subset. NOTE: assumes not < and not > --> ==
    auto lb = std::lower_bound(first, last, key);
    return (lb == end() || *lb != key) ? end() : lb;
  }
  const_iterator overlap_lower_bound(
      const ParentT& parent, position_t coordinate) const {
    // offsets into elements_ that are first vs last for given parent
    size_t idx_first = parent_idx_offsets_[parent.idx_];
    size_t idx_last = parent_idx_offsets_[1 + parent.idx_];
    // first interval that can overlap with coordinate has end >= coordinate
    if constexpr(HAS_OVERLAPS) {
      auto cummax_lb = std::lower_bound(
          elements_end_cummax_.begin() + idx_first,
          elements_end_cummax_.begin() + idx_last,
          coordinate);
      return begin() + (cummax_lb - elements_end_cummax_.begin());
    } else {
      return std::lower_bound(
          begin() + idx_first,
          begin() + idx_last,
          coordinate,
          [](const RegionT& region, const position_t& x) {
          return region.coordinates.last_pos() < x;
          });
    }
  }
  const_iterator overlap_upper_bound(
      const ParentT& parent, position_t coordinate) const {
    // first interval that can overlap with coordinate has start > coordinate
    return std::upper_bound(
        // index into elements that share parent
        begin() + parent_idx_offsets_[parent.idx_],
        begin() + parent_idx_offsets_[1 + parent.idx_],
        coordinate,
        [](const RegionT& region, const position_t& x) {
        return region.coordinates.first_pos() < x;
        });
  }

  explicit Regions(vecT&& x)
      : elements_{std::move(x)},
        parents_{elements_.empty() ? nullptr : elements_[0].parent().ptr_},
        parent_idx_offsets_{validate_and_parent_offsets(elements_, parents_)},
        elements_end_cummax_{get_cummax_ends(elements_)} { }
  Regions() : Regions{vecT{}} { }
  Regions(const Regions&) = default;
  Regions(Regions&&) = default;
  Regions& operator=(const Regions&) = delete;
  Regions& operator=(Regions&&) = delete;

  friend inline bool operator==(const Regions& x, const Regions& y) {
    return std::tie(x.parents_, x.elements_)
      == std::tie(y.parents_, y.elements_);
  }
};

}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_REGIONS_HPP
