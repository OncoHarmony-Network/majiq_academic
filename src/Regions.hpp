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

// compare regions with respect to contig-coordinates, including strand
// NOTE: Satisfies Compare (full ordering) when between same type
template <class T1, class T2>
struct CompareContigStranded {
  inline bool operator()(
      const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
    const Gene& gx = x.gene.get();
    const Gene& gy = y.gene.get();
    // still complete ordering, just gene least important
    return std::tie(gx.contig, x.coordinates, gx.strand, x.gene)
      < std::tie(gy.contig, y.coordinates, gy.strand, y.gene);
  }
  inline bool operator()(
      const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
    return x < y;
  }
  inline bool operator()(
      const GeneRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
    return x < y;
  }
  inline bool operator()(
      const ContigRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
    return x < y;
  }
};

// unstranded comparison with contig-coordinates (ignore strand, gene)
// NOTE: partial ordering (not Compare requirements)
// Compatible with CompareContigStranded -- sorted by CompareContigStranded -->
// sorted by CompareContigUnstranded
template <class T1, class T2>
struct CompareContigUnstranded {
  inline bool operator()(
      const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
    const Gene& gx = x.gene.get();
    const Gene& gy = y.gene.get();
    return std::tie(gx.contig, x.coordinates)
      < std::tie(gy.contig, y.coordinates);
  }
  inline bool operator()(
      const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
    return std::tie(x.contig, x.coordinates)
      < std::tie(y.contig, y.coordinates);
  }
  inline bool operator()(
      const ContigRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
    const Gene& gy = y.gene.get();
    return std::tie(x.contig, x.coordinates)
      < std::tie(gy.contig, y.coordinates);
  }
  inline bool operator()(
      const GeneRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
    const Gene& gx = x.gene.get();
    return std::tie(gx.contig, x.coordinates)
      < std::tie(y.contig, y.coordinates);
  }
};

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
template <class RegionT, class CompareRegionT = std::less<RegionT>>
class GeneRegions {
 public:
  using value_type = RegionT;
  using key_compare = CompareRegionT;
  using vecRegionT = std::vector<RegionT>;
  using _setRegionT_unordered = std::unordered_set<RegionT>;
  using _setRegionT_ordered = std::set<RegionT, CompareRegionT>;
  using setRegionT = _setRegionT_unordered;  // what we are actually using
  using IntervalT = typename RegionT::IntervalT;
  static_assert(std::is_base_of<GeneRegion<IntervalT>, RegionT>::value,
      "GeneRegions type must be subclass of GeneRegion");

 protected:
  vecRegionT sorted_vec_;
  setRegionT region_set_;  // temporary regions that couldn't be added/sorted

  // static functions used by constructors, etc.
  /**
   * Combine two containers of unique RegionT into single sorted unique vector
   *
   * @note container other assumed to be sorted if matches _setRegionT_ordered
   * or other_sorted is True
   */
  template <bool other_sorted = false, class T>
  static vecRegionT _combine_sorted(vecRegionT sorted, T other) {
    static_assert(
        std::is_same<RegionT, typename T::value_type>::value,
        "other must have RegionT values in _combine_sorted");
    if (other.size() == 0) {
      return sorted;
    }
    vecRegionT result;
    result.reserve(sorted.size() + other.size());
    if constexpr(other_sorted || std::is_same<T, _setRegionT_ordered>::value) {
      // region_set ordered the way we need it to be already
      std::set_union(
          sorted.begin(), sorted.end(), other.begin(), other.end(),
          std::back_inserter(result), CompareRegionT{});
    } else {
      // not sure if sorted, so dump to temporary array to sort the way we want
      vecRegionT other_vec{other.begin(), other.end()};
      std::sort(other_vec.begin(), other_vec.end(), CompareRegionT{});
      std::set_union(
          sorted.begin(), sorted.end(), other_vec.begin(), other_vec.end(),
          std::back_inserter(result), CompareRegionT{});
    }
    return result;
  }

  /**
   * get sorted vector, remapping genes in the regions passed in from two
   * containers
   */
  template <class T1, class T2>
  static vecRegionT _remap_combine(
      T1 regions1, T2 regions2,
      const std::shared_ptr<Genes>& new_genes) {
    static_assert(
        std::is_same<RegionT, typename T1::value_type>::value,
        "regions1 must have elements of RegionT in _remap_combine");
    static_assert(
        std::is_same<RegionT, typename T2::value_type>::value,
        "regions2 must have elements of RegionT in _remap_combine");
    // initialize/reserve resulting vector
    vecRegionT out;
    out.reserve(regions1.size() + regions2.size());
    // transform exons using genes, add to back of new_vec
    auto remapper = [&new_genes](RegionT x) -> RegionT {
      x.gene = x.gene.remapped(new_genes);
      return x;
    };
    std::transform(
        regions1.begin(), regions1.end(), std::back_inserter(out), remapper);
    std::transform(
        regions2.begin(), regions2.end(), std::back_inserter(out), remapper);
    // sort new_vec
    std::sort(out.begin(), out.end(), CompareRegionT{});
    // return the result
    return out;
  }

  // helper functions for determining if a region is present
  bool _sorted_contains(const RegionT& x) const {
    auto eq_range = std::equal_range(
        sorted_vec_.begin(), sorted_vec_.end(), x, CompareRegionT{});
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
    if (region_set_.size() > 0) {
      sorted_vec_ = _remap_combine(sorted_vec_, region_set_, new_genes);
      region_set_.clear();
    } else {
      // remap/sort in place...
      std::for_each(sorted_vec_.begin(), sorted_vec_.end(),
          [&new_genes](RegionT& x) { x.gene = x.gene.remapped(new_genes); });
      std::sort(sorted_vec_.begin(), sorted_vec_.end(), CompareRegionT{});
    }
    return;
  }
  // access exons
  const RegionT& operator[](size_t idx) {
    sort();
    return sorted_vec_[idx];
  }
  const RegionT& operator[](size_t idx) const {
    if (unsorted()) {
      std::cerr << "operator[] on unsorted const Regions not intended/ideal\n";
      return GeneRegions<RegionT, CompareRegionT>(this)[idx];
    } else {
      return sorted_vec_[idx];
    }
  }
  // iterators
  constexpr typename vecRegionT::const_iterator begin() {
    sort();
    return sorted_vec_.begin();
  }
  constexpr typename vecRegionT::const_iterator begin() const {
    if (unsorted()) {
      std::cerr << "begin() on unsorted const Regions not intended/ideal\n";
      return GeneRegions<RegionT, CompareRegionT>(this).begin();
    } else {
      return sorted_vec_.begin();
    }
  }
  constexpr typename vecRegionT::const_iterator end() {
    sort();
    return sorted_vec_.end();
  }
  constexpr typename vecRegionT::const_iterator end() const {
    if (unsorted()) {
      std::cerr << "end() on unsorted const Regions not intended/ideal\n";
      return GeneRegions<RegionT, CompareRegionT>(this).end();
    } else {
      return sorted_vec_.end();
    }
  }

  // membership
  bool contains(const RegionT& x) const {
    return _sorted_contains(x) || _set_contains(x);
  }

  // insertion of new elements
  void insert(const RegionT& x) {
    if (sorted_vec_.size() == 0 || CompareRegionT{}(sorted_vec_.back(), x)) {
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
