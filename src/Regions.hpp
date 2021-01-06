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

#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
namespace detail {

template <class T>
struct ContigRegion {
  using IntervalT = T;
  static_assert(std::is_base_of<Interval, IntervalT>::value,
      "IntervalT must be derived from Interval (Open or Closed)");
 public:
  // location
  KnownContig contig;
  IntervalT coordinates;
  GeneStrandness strand;

  // constructors
  ContigRegion(KnownContig _contig, IntervalT _coordinates,
      GeneStrandness _strand)
      : contig{_contig}, coordinates{_coordinates}, strand{_strand} {
  }
  ContigRegion(const ContigRegion& x)
      : ContigRegion{x.contig, x.coordinates, x.strand} {
  }
};

template <class T>
struct GeneRegion {
  using IntervalT = T;
  static_assert(std::is_base_of<Interval, IntervalT>::value,
      "IntervalT must be derived from Interval (Open or Closed)");
 public:
  // location
  KnownGene gene;
  IntervalT coordinates;

  // constructors
  GeneRegion(KnownGene _gene, IntervalT _coordinates)
      : gene{_gene}, coordinates{_coordinates} {
  }
  GeneRegion(const GeneRegion& x) : GeneRegion{x.gene, x.coordinates} {}
};

// order regions by genomic position and strand
template <class T1, class T2>
inline bool operator<(
    const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
/**
 * Order regions by gene, then position
 */
template <class T1, class T2>
inline bool operator<(
    const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return std::tie(x.gene, x.coordinates) < std::tie(y.gene, y.coordinates);
}
/**
 * Compare to genes
 */
template <class T>
inline bool operator<(const GeneRegion<T>& lhs, const KnownGene& rhs) noexcept {
  return lhs.gene < rhs;
}
template <class T>
inline bool operator<(const KnownGene& lhs, const GeneRegion<T>& rhs) noexcept {
  return lhs < rhs.gene;
}
/**
 * Compare to contig regions
 */
template <class T1, class T2>
inline bool operator<(
    const GeneRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  const Gene& gx = x.gene.get();
  return std::tie(gx.contig, x.coordinates, gx.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
template <class T1, class T2>
inline bool operator<(
    const ContigRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  const Gene& gy = y.gene.get();
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(gy.contig, y.coordinates, gy.strand);
}

// compare regions with respect to contig-coordinates but also stranded
// (when involving ContigRegion just an alias for operator<)
template <class T1, class T2>
inline bool CompareContigStranded(
    const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  const Gene& gx = x.gene.get();
  const Gene& gy = y.gene.get();
  return std::tie(gx.contig, x.coordinates, gx.strand)
    < std::tie(gy.contig, y.coordinates, gy.strand);
}
template <class T1, class T2>
inline bool CompareContigStranded(
    const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return x < y;
}
template <class T1, class T2>
inline bool CompareContigStranded(
    const GeneRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return x < y;
}
template <class T1, class T2>
inline bool CompareContigStranded(
    const ContigRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return x < y;
}

// unstranded comparison with contig-coordinates
template <class T1, class T2>
inline bool CompareContigUnstranded(
    const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  const Gene& gx = x.gene.get();
  const Gene& gy = y.gene.get();
  return std::tie(gx.contig, x.coordinates)
    < std::tie(gy.contig, y.coordinates);
}
template <class T1, class T2>
inline bool CompareContigUnstranded(
    const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return std::tie(x.contig, x.coordinates) < std::tie(y.contig, y.coordinates);
}
template <class T1, class T2>
inline bool CompareContigUnstranded(
    const ContigRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  const Gene& gy = y.gene.get();
  return std::tie(x.contig, x.coordinates) < std::tie(gy.contig, y.coordinates);
}
template <class T1, class T2>
inline bool CompareContigUnstranded(
    const GeneRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  const Gene& gx = x.gene.get();
  return std::tie(gx.contig, x.coordinates) < std::tie(y.contig, y.coordinates);
}

// allow regions to be passed into output stream (e.g. std::cout)
template <class T>
std::ostream& operator<<(std::ostream& os, const ContigRegion<T>& x) noexcept {
  os << x.contig.get()
    << ":" << x.strand
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
}
template <class T>
std::ostream& operator<<(std::ostream& os, const GeneRegion<T>& x) noexcept {
  os << x.gene.get()
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
}

// derived comparisons (ContigRegion, ContigRegion)
template <class T1, class T2>
inline bool operator>(
    const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return y < x;
}
template <class T1, class T2>
inline bool operator<=(
    const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return !(y < x);
}
template <class T1, class T2>
inline bool operator>=(
    const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return !(x < y);
}

// derived comparisons (GeneRegion, GeneRegion)
template <class T1, class T2>
inline bool operator>(
    const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return y < x;
}
template <class T1, class T2>
inline bool operator>=(
    const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return !(x < y);
}
template <class T1, class T2>
inline bool operator<=(
    const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return !(y < x);
}

// derived comparisons (KnownGene, GeneRegion)
template <class T>
inline bool operator>(const GeneRegion<T>& x, const KnownGene& y) noexcept {
  return y < x;
}
template <class T>
inline bool operator>(const KnownGene& x, const GeneRegion<T>& y) noexcept {
  return y < x;
}
template <class T>
inline bool operator<=(const GeneRegion<T>& x, const KnownGene& y) noexcept {
  return !(y < x);
}
template <class T>
inline bool operator<=(const KnownGene& x, const GeneRegion<T>& y) noexcept {
  return !(y < x);
}
template <class T>
inline bool operator>=(const GeneRegion<T>& x, const KnownGene& y) noexcept {
  return !(x < y);
}
template <class T>
inline bool operator>=(const KnownGene& x, const GeneRegion<T>& y) noexcept {
  return !(x < y);
}

// derived comparisons (ContigRegion, GeneRegion)
template <class T1, class T2>
inline bool operator>(
    const GeneRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return y < x;
}
template <class T1, class T2>
inline bool operator>(
    const ContigRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return y < x;
}
template <class T1, class T2>
inline bool operator<=(
    const GeneRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return !(y < x);
}
template <class T1, class T2>
inline bool operator<=(
    const ContigRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return !(y < x);
}
template <class T1, class T2>
inline bool operator>=(
    const GeneRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return !(x < y);
}
template <class T1, class T2>
inline bool operator>=(
    const ContigRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return !(x < y);
}


// implement template container of gene regions
template <class RegionT>
class GeneRegions {
 public:
  using vecRegionT = std::vector<RegionT>;
  using setRegionT = std::unordered_set<RegionT>;
  using IntervalT = typename RegionT::IntervalT;
  static_assert(std::is_base_of<GeneRegion<IntervalT>, RegionT>::value,
      "GeneRegions type must be subclass of GeneRegion");

 protected:
  vecRegionT sorted_vec_;
  setRegionT region_set_;  // temporary regions that couldn't be added/sorted

  // static functions used by constructors, etc.
  /**
   * Get sorted vector out of sorted vector and new set of items
   * @note explicitly defining for unordered_set vs (ordered) set, which
   * doesn't require sorting step
   */
  static vecRegionT _combine_sorted(
      vecRegionT sorted, std::unordered_set<RegionT> unsorted) {
    // get unsorted values in sorted order
    vecRegionT unsorted_vec{unsorted.begin(), unsorted.end()};
    std::sort(unsorted_vec.begin(), unsorted_vec.end());
    // merge the two sorted vectors for final result
    vecRegionT result;
    result.reserve(sorted.size() + unsorted.size());
    std::merge(
        sorted.begin(), sorted.end(),
        unsorted_vec.begin(), unsorted_vec.end(),
        std::back_inserter(result));
    // return final result
    return result;
  }
  /**
   * Get sorted vector out of sorted vector and sorted set of items
   * @note explicitly defining for unordered_set vs (ordered) set, which
   * doesn't require sorting step
   */
  static vecRegionT _combine_sorted(
      vecRegionT sorted, std::set<RegionT> set_sorted) {
    vecRegionT result;
    result.reserve(sorted.size() + set_sorted.size());
    std::merge(
        sorted.begin(), sorted.end(),
        set_sorted.begin(), set_sorted.end(),
        std::back_inserter(result));
    // return final result
    return result;
  }

  /**
   * sorted vector when remapping genes in the exons
   */
  static vecRegionT _remap_combine(
      vecRegionT regions_vec, setRegionT regions_set,
      const std::shared_ptr<Genes>& new_genes) {
    // initialize/reserve resulting vector
    vecRegionT new_vec;
    new_vec.reserve(regions_vec.size() + regions_set.size());
    // transform exons using genes, add to back of new_vec
    auto remapper = [new_genes](RegionT x) -> RegionT {
      x.gene = x.gene.remapped(new_genes);
      return x;
    };
    std::transform(
        regions_vec.begin(), regions_vec.end(),
        std::back_inserter(new_vec), remapper);
    std::transform(
        regions_set.begin(), regions_set.end(),
        std::back_inserter(new_vec), remapper);
    // sort new_vec
    std::sort(new_vec.begin(), new_vec.end());
    // return the result
    return new_vec;
  }

  // helper functions for determining if a region is present
  bool _sorted_contains(const RegionT& x) const {
    auto eq_range = std::equal_range(sorted_vec_.begin(), sorted_vec_.end(), x);
    return eq_range.first != eq_range.second;
  }
  bool _set_contains(const RegionT& x) const {
    return region_set_.count(x) > 0;
  }

 public:
  size_t size() const { return sorted_vec_.size() + region_set_.size(); }
  bool unsorted() const { return region_set_.size() > 0; }
  void sort() {
    if (unsorted()) {
      sorted_vec_ = _combine_sorted(sorted_vec_, region_set_);
      region_set_.clear();
    }
    return;
  }
  void remap_genes(const std::shared_ptr<Genes>& new_genes) {
    sorted_vec_ = _remap_combine(sorted_vec_, region_set_, new_genes);
    region_set_.clear();
    return;
  }
  // access exons
  const RegionT& operator[](size_t idx) {
    sort();
    return sorted_vec_[idx];
  }
  const RegionT& operator[](size_t idx) const {
    // NOTE assumes sorted, so make sure constructors always sorted
    return sorted_vec_[idx];
  }

  // membership
  bool contains(const RegionT& x) const {
    return _sorted_contains(x) || _set_contains(x);
  }

  // insertion of new elements
  void insert(const RegionT& x) {
    if (sorted_vec_.size() == 0 || sorted_vec_.back() < x) {
      // it belongs at end of the sorted vector
      sorted_vec_.push_back(x);
      return;
    } else if (contains(x)) {
      // we already have it
      return;
    } else {
      // can't put it at end of sorted vector, but we don't have it
      region_set_.insert(x);
      return;
    }
  }

  // constructors
  GeneRegions() : sorted_vec_{}, region_set_{} { }
  GeneRegions(const vecRegionT& v, const setRegionT& s)
      : sorted_vec_{_combine_sorted(v, s)}, region_set_{} {
  }
  explicit GeneRegions(const vecRegionT& v) : GeneRegions{v, setRegionT{}} { }
  GeneRegions(const GeneRegions& x)
      : GeneRegions{x.sorted_vec_, x.region_set_} {
  }
  GeneRegions(const vecRegionT& v, const setRegionT& s,
      const std::shared_ptr<Genes>& g)
      : sorted_vec_{_remap_combine(v, s, g)}, region_set_{} {
  }
  GeneRegions(const vecRegionT& v, const std::shared_ptr<Genes>& g)
      : GeneRegions{v, setRegionT{}, g} {
  }
  GeneRegions(const GeneRegions& x, const std::shared_ptr<Genes>& g)
      : GeneRegions{x.sorted_vec_, x.region_set_, g} {
  }

  // comparisons
  friend inline bool operator==(const GeneRegions& x, const GeneRegions& y) {
    return x.sorted_vec_ == y.sorted_vec_ && x.region_set_ == y.region_set_;
  }

  // view into regionT vector for pybind11
  const vecRegionT& data() { sort(); return sorted_vec_; }
};
}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_REGIONS_HPP
