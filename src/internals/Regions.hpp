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


template <
  class RegionT,
  class CompareRegionT = std::less<RegionT>>
class GeneRegions {
 public:
  using IntervalT = typename RegionT::IntervalT;
  using DataT = typename RegionT::DataT;
  using BaseRegion = GeneRegion<IntervalT, DataT>;
  static_assert(std::is_base_of<BaseRegion, RegionT>::value,
      "GeneRegions type must be subclass of GeneRegion");
  using value_type = RegionT;
  using key_compare = CompareRegionT;

  using vecT = std::vector<value_type>;
  using vecConstT = std::vector<const value_type>;
  using const_iterator = typename vecT::const_iterator;
  struct NoCheckValid { };

 private:
  vecT elements_;
  std::shared_ptr<Genes> genes_;

  static bool is_valid(const vecT& x) {
    return x.end() == std::adjacent_find(x.begin(), x.end(),
        [](const value_type& a, const value_type& b) {
        // should be sorted and have equal pointers, return if not
        return !(
            CompareRegionT{}(b, a)
            || static_cast<BaseRegion>(a) == static_cast<BaseRegion>(b)
            || a.gene.known_genes == b.gene.known_genes);
        });
  }

 public:
  GeneRegions() : elements_{}, genes_{nullptr} { }
  GeneRegions(vecT&& x, NoCheckValid)
      : elements_{x}, genes_{x.empty() ? nullptr : x[0].gene.known_genes} { }
  GeneRegions(const vecT& x, NoCheckValid)
      : elements_{x}, genes_{x.empty() ? nullptr : x[0].gene.known_genes} { }
  explicit GeneRegions(vecT&& x) : GeneRegions{x, NoCheckValid{}} {
    if (!is_valid(x)) {
      throw std::invalid_argument("vector input to GeneRegions is invalid");
    }
  }
  explicit GeneRegions(const vecT& x) : GeneRegions{x, NoCheckValid{}} {
    if (!is_valid(x)) {
      throw std::invalid_argument("vector input to GeneRegions is invalid");
    }
  }
  GeneRegions(const GeneRegions& x) = default;
  GeneRegions(GeneRegions&& x) = default;
  GeneRegions& operator=(const GeneRegions& x) = default;
  GeneRegions& operator=(GeneRegions&& x) = default;

  size_t size() const noexcept { return elements_.size(); }
  bool empty() const noexcept { return elements_.empty(); }

  const value_type& operator[](size_t idx) const { return elements_[idx]; }
  const_iterator begin() const { return elements_.cbegin(); }
  const_iterator end() const { return elements_.cend(); }
  const vecT& data() { return elements_; }
  const vecConstT& data() const { return elements_; }
  const std::shared_ptr<Genes>& genes() const { return genes_; }

  const_iterator find(const value_type& key) const {
    // run lower-bound and uppper-bound simultaneously
    auto&& [lb, ub] = std::equal_range(begin(), end(), key, key_compare{});
    // end of range if no match (lb==ub), otherwise first match
    return lb == ub ? end() : lb;
  }
  size_t count(const value_type& key) const {
    return find(key) == end() ? 0 : 1;
  }
  friend inline bool operator==(const GeneRegions& x, const GeneRegions& y) {
    return std::tie(x.genes_, x.elements_) == std::tie(y.genes_, y.elements_);
  }
};

template <class RegionT>
class ContigRegions {
 public:
  using IntervalT = typename RegionT::IntervalT;
  using DataT = typename RegionT::DataT;
  using BaseRegion = ContigRegion<IntervalT, DataT>;
  static_assert(std::is_base_of<BaseRegion, RegionT>::value,
      "ContigRegions type must be subclass of ContigRegion");
  using value_type = RegionT;

  using vecT = std::vector<value_type>;
  using vecConstT = std::vector<const value_type>;
  using const_iterator = typename vecT::const_iterator;
  struct NoCheckValid { };

 protected:
  vecT elements_;
  std::shared_ptr<Contigs> contigs_;

  static bool is_valid(const vecT& x) {
    return x.end() == std::adjacent_find(x.begin(), x.end(),
        [](const value_type& a, const value_type& b) {
        // should be sorted and have equal pointers, return if not
        return !(a < b && a.contig.ptr_ == b.contig.ptr_);
        });
  }

 public:
  ContigRegions() : elements_{}, contigs_{nullptr} { }
  ContigRegions(vecT&& x, NoCheckValid)
      : elements_{x},
        contigs_{x.empty() ? nullptr : x[0].contig.ptr_} { }
  ContigRegions(const vecT& x, NoCheckValid)
      : elements_{x},
        contigs_{x.empty() ? nullptr : x[0].contig.ptr_} { }
  explicit ContigRegions(vecT&& x) : ContigRegions{x, NoCheckValid{}} {
    if (!is_valid(x)) {
      throw std::invalid_argument("vector input to ContigRegions is invalid");
    }
  }
  explicit ContigRegions(const vecT& x) : ContigRegions{x, NoCheckValid{}} {
    if (!is_valid(x)) {
      throw std::invalid_argument("vector input to ContigRegions is invalid");
    }
  }
  ContigRegions(const ContigRegions& x) = default;
  ContigRegions(ContigRegions&& x) = default;
  ContigRegions& operator=(const ContigRegions& x) = default;
  ContigRegions& operator=(ContigRegions&& x) = default;

  size_t size() const noexcept { return elements_.size(); }
  bool empty() const noexcept { return elements_.empty(); }

  const value_type& operator[](size_t idx) const { return elements_[idx]; }
  const_iterator begin() const { return elements_.cbegin(); }
  const_iterator end() const { return elements_.cend(); }
  const vecT& data() { return elements_; }
  const vecConstT& data() const { return elements_; }
  const std::shared_ptr<Contigs>& contigs() const { return contigs_; }

  template<class Compare>
  const_iterator find(const value_type& key) const {
    // run lower-bound and uppper-bound simultaneously
    auto&& [lb, ub] = std::equal_range(begin(), end(), key, Compare{});
    // end of range if no match (lb==ub), otherwise first match
    return lb == ub ? end() : lb;
  }
  size_t count(const value_type& key) const {
    return find(key) == end() ? 0 : 1;
  }
  friend inline bool operator==(
      const ContigRegions& x, const ContigRegions& y) {
    return std::tie(x.contigs_, x.elements_)
      == std::tie(y.contigs_, y.elements_);
  }
};
}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_REGIONS_HPP
