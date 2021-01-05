/**
 * Exons.hpp
 *
 * Exons for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_EXONS_HPP
#define MAJIQ_EXONS_HPP

#include <tuple>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <memory>
#include <boost/functional/hash.hpp>

#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
struct Exon : detail::GeneRegion {
 public:
  // constructors
  Exon(KnownGene _gene, ClosedInterval _coordinates)
      : detail::GeneRegion{_gene, _coordinates} {
  }
  Exon(const Exon& x) : Exon{x.gene, x.coordinates} {}

  // comparison for equality
  bool operator==(const Exon& rhs) const {
    return std::tie(gene, coordinates)
      == std::tie(rhs.gene, rhs.coordinates);
  }
};

// override boost hashing
std::size_t hash_value(const Exon& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// specialize std::hash for Exon
namespace std {
template <> struct hash<majiq::Exon> {
  std::size_t operator()(const majiq::Exon& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
/**
 * Exons in sorted contiguous container, with extra hash table for exons that
 * were added in unsorted order to defer sorting until necessary
 */
class Exons {
 public:
  using exon_set_t = std::unordered_set<Exon>;

 private:
  // sorted vector of exons
  std::vector<Exon> exons_vec_;
  // unsorted members
  exon_set_t exons_set_;

  // get a sorted vector out of a sorted vector and unordered set
  static std::vector<Exon> _combine_sorted(
      std::vector<Exon> sorted, std::unordered_set<Exon> unsorted) {
    // get unsorted values in sorted order
    std::vector<Exon> unsorted_vec{unsorted.begin(), unsorted.end()};
    std::sort(unsorted_vec.begin(), unsorted_vec.end());
    // merge the two sorted vectors for final result
    std::vector<Exon> result;
    result.reserve(sorted.size() + unsorted.size());
    std::merge(
        sorted.begin(), sorted.end(),
        unsorted_vec.begin(), unsorted_vec.end(),
        std::back_inserter(result));
    // return final result
    return result;
  }
  static std::vector<Exon> _combine_sorted(
      std::vector<Exon> sorted, std::set<Exon> set_sorted) {
    std::vector<Exon> result;
    result.reserve(sorted.size() + set_sorted.size());
    std::merge(
        sorted.begin(), sorted.end(),
        set_sorted.begin(), set_sorted.end(),
        std::back_inserter(result));
    // return final result
    return result;
  }
  // get sorted vector when remapping genes in the exons
  static std::vector<Exon> _remap_combine(
      std::vector<Exon> exon_vec, exon_set_t exon_set,
      const std::shared_ptr<Genes>& new_genes) {
    // initialize/reserve resulting vector
    std::vector<Exon> new_vec;
    new_vec.reserve(exon_vec.size() + exon_set.size());
    // transform exons using genes, add to back of new_vec
    auto remapper = [new_genes](Exon x) -> Exon {
      x.gene = x.gene.remapped(new_genes);
      return x;
    };
    std::transform(
        exon_vec.begin(), exon_vec.end(),
        std::back_inserter(new_vec), remapper);
    std::transform(
        exon_set.begin(), exon_set.end(),
        std::back_inserter(new_vec), remapper);
    // sort new_vec
    std::sort(new_vec.begin(), new_vec.end());
    // return the result
    return new_vec;
  }

  bool _sorted_contains(const Exon& x) const {
    auto eq_range = std::equal_range(exons_vec_.begin(), exons_vec_.end(), x);
    return eq_range.first != eq_range.second;
  }
  bool _unsorted_contains(const Exon& x) const {
    return exons_set_.count(x) > 0;
  }

 public:
  // how many elements we have
  size_t size() const { return exons_vec_.size() + exons_set_.size(); }
  // manage if it's sorted
  bool unsorted() const { return exons_set_.size() > 0; }
  void sort() {
    if (unsorted()) {
      exons_vec_ = _combine_sorted(exons_vec_, exons_set_);
      exons_set_.clear();  // all values sorted now
    }
  }
  // remap genes and sort
  void remap_genes(const std::shared_ptr<Genes>& new_genes) {
    exons_vec_ = _remap_combine(exons_vec_, exons_set_, new_genes);
    exons_set_.clear();  // all values sorted now
  }
  // access exons
  const Exon& operator[](size_t idx) {
    sort();
    return exons_vec_[idx];
  }
  /*
   * @note assumes container is sorted -- cannot require if const
   */
  const Exon& operator[](size_t idx) const { return exons_vec_[idx]; }

  // membership, insertion
  bool contains(const Exon& x) const {
    return _sorted_contains(x) || _unsorted_contains(x);
  }
  void insert(const Exon& x) {
    if (exons_vec_.size() == 0 || exons_vec_.back() < x) {
      // it belongs at end of the sorted vector
      exons_vec_.push_back(x);
      return;
    } else if (contains(x)) {
      // we already have it
      return;
    } else {
      // can't put it at end of sorted vector, but we don't have it
      exons_set_.insert(x);
      return;
    }
  }

  // constructors
  Exons() : exons_vec_{}, exons_set_{} { }
  Exons(const std::vector<Exon>& v, const exon_set_t& s)
      : exons_vec_{_combine_sorted(v, s)}, exons_set_{} {
  }
  explicit Exons(const std::vector<Exon>& v) : Exons{v, exon_set_t{}} { }
  Exons(const Exons& x) : Exons{x.exons_vec_, x.exons_set_} { }
  Exons(const std::vector<Exon>& v, const exon_set_t& s,
      std::shared_ptr<Genes> g)
      : exons_vec_{_remap_combine(v, s, g)}, exons_set_{} {
  }
  Exons(const std::vector<Exon>& v, std::shared_ptr<Genes> g)
      : Exons{v, exon_set_t{}, g} {
  }
  Exons(const Exons& x, std::shared_ptr<Genes> g)
      : Exons{x.exons_vec_, x.exons_set_, g} {
  }

  // comparisons
  friend inline bool operator==(const Exons& x, const Exons& y);
};

inline bool operator==(const Exons& x, const Exons& y) {
  return x.exons_vec_ == y.exons_vec_ && x.exons_set_ == y.exons_set_;
}
}  // namespace majiq

#endif  // MAJIQ_EXONS_HPP
