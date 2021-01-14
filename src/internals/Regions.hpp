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

#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"


namespace majiq {
namespace detail {

// empty data type for carrying data associated with regions
struct EmptyDataT {
  EmptyDataT() = default;
  EmptyDataT(const EmptyDataT& x) = default;
  EmptyDataT(EmptyDataT&& x) = default;
  EmptyDataT& operator=(const EmptyDataT& x) = default;
  EmptyDataT& operator=(EmptyDataT&& x) = default;
};

template <class T, class DT = EmptyDataT>
struct ContigRegion {
 public:
  using IntervalT = T;
  using DataT = DT;
  static_assert(std::is_base_of<Interval, IntervalT>::value,
      "IntervalT must be derived from Interval (Open or Closed)");

  // location
  KnownContig contig;
  IntervalT coordinates;
  GeneStrandness strand;
  DataT data;

  // constructors
  ContigRegion(
      KnownContig _contig, IntervalT _coordinates, GeneStrandness _strand,
      DataT _data)
      : contig{_contig},
        coordinates{_coordinates},
        strand{_strand},
        data{_data} {
  }
  ContigRegion(KnownContig _contig, IntervalT _coordinates,
      GeneStrandness _strand)
      : ContigRegion{_contig, _coordinates, _strand, DataT{}} {
  }
  ContigRegion()
      : ContigRegion{KnownContig{}, IntervalT{}, GeneStrandness::AMBIGUOUS,
        DataT{}} {
  }
  ContigRegion(const ContigRegion& x) = default;
  ContigRegion(ContigRegion&& x) = default;
  ContigRegion& operator=(const ContigRegion& x) = default;
  ContigRegion& operator=(ContigRegion&& x) = default;

  // drop data between templates
  template <class OtherDataT>
  ContigRegion(const ContigRegion<IntervalT, OtherDataT>& x)
      : ContigRegion{x.contig, x.coordinates, x.strand} { }
  template <class OtherDataT>
  ContigRegion(ContigRegion<IntervalT, OtherDataT>&& x)
      : ContigRegion{x.contig, x.coordinates, x.strand} { }
};

template <class T, class DT = EmptyDataT>
struct GeneRegion {
  using IntervalT = T;
  using DataT = DT;
  static_assert(std::is_base_of<Interval, IntervalT>::value,
      "IntervalT must be derived from Interval (Open or Closed)");
 public:
  // location
  KnownGene gene;
  IntervalT coordinates;
  DataT data;

  // constructors
  GeneRegion(KnownGene _gene, IntervalT _coordinates, DataT _data)
      : gene{_gene}, coordinates{_coordinates}, data{_data} {
  }
  GeneRegion(KnownGene _gene, IntervalT _coordinates)
      : GeneRegion{_gene, _coordinates, DataT{}} {
  }
  GeneRegion() : GeneRegion{KnownGene{}, IntervalT{}} { }
  GeneRegion(const GeneRegion& x) = default;
  GeneRegion(GeneRegion&& x) = default;
  GeneRegion& operator=(const GeneRegion& x) = default;
  GeneRegion& operator=(GeneRegion&& x) = default;
};

// order regions by genomic position and strand
template <class T1, class T2, class D1, class D2>
inline bool operator<(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
template <class T1, class T2, class D1, class D2>
inline bool operator==(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    == std::tie(y.contig, y.coordinates, y.strand);
}
/**
 * Order regions by gene, then position
 */
template <class T1, class T2, class D1, class D2>
inline bool operator<(
    const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return std::tie(x.gene, x.coordinates) < std::tie(y.gene, y.coordinates);
}
template <class T1, class T2, class D1, class D2>
inline bool operator==(
    const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return std::tie(x.gene, x.coordinates) == std::tie(y.gene, y.coordinates);
}
/**
 * Compare to genes
 */
template <class T, class D>
inline bool operator<(
    const GeneRegion<T, D>& lhs, const KnownGene& rhs) noexcept {
  return lhs.gene < rhs;
}
template <class T, class D>
inline bool operator<(
    const KnownGene& lhs, const GeneRegion<T, D>& rhs) noexcept {
  return lhs < rhs.gene;
}
/**
 * Compare to contig regions
 */
template <class T1, class T2, class D1, class D2>
inline bool operator<(
    const GeneRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  const Gene& gx = x.gene.get();
  return std::tie(gx.contig, x.coordinates, gx.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
template <class T1, class T2, class D1, class D2>
inline bool operator<(
    const ContigRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  const Gene& gy = y.gene.get();
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(gy.contig, y.coordinates, gy.strand);
}

// compare regions with respect to contig-coordinates, including strand
// NOTE: Satisfies Compare (full ordering) when between same type
template <class T1, class T2, class D1, class D2>
struct CompareContigStranded {
  inline bool operator()(
      const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
    const Gene& gx = x.gene.get();
    const Gene& gy = y.gene.get();
    // still complete ordering, just gene least important
    return std::tie(gx.contig, x.coordinates, gx.strand, x.gene)
      < std::tie(gy.contig, y.coordinates, gy.strand, y.gene);
  }
  inline bool operator()(
      const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
    return x < y;
  }
  inline bool operator()(
      const GeneRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
    return x < y;
  }
  inline bool operator()(
      const ContigRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
    return x < y;
  }
};

// unstranded comparison with contig-coordinates (ignore strand, gene)
// NOTE: partial ordering (not Compare requirements)
// Compatible with CompareContigStranded -- sorted by CompareContigStranded -->
// sorted by CompareContigUnstranded
template <class T1, class T2, class D1, class D2>
struct CompareContigUnstranded {
  inline bool operator()(
      const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
    const Gene& gx = x.gene.get();
    const Gene& gy = y.gene.get();
    return std::tie(gx.contig, x.coordinates)
      < std::tie(gy.contig, y.coordinates);
  }
  inline bool operator()(
      const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
    return std::tie(x.contig, x.coordinates)
      < std::tie(y.contig, y.coordinates);
  }
  inline bool operator()(
      const ContigRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
    const Gene& gy = y.gene.get();
    return std::tie(x.contig, x.coordinates)
      < std::tie(gy.contig, y.coordinates);
  }
  inline bool operator()(
      const GeneRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
    const Gene& gx = x.gene.get();
    return std::tie(gx.contig, x.coordinates)
      < std::tie(y.contig, y.coordinates);
  }
};

// allow regions to be passed into output stream (e.g. std::cout)
template <class T, class D>
std::ostream& operator<<(
    std::ostream& os, const ContigRegion<T, D>& x) noexcept {
  os << x.contig.get()
    << ":" << x.strand
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
}
template <class T, class D>
std::ostream& operator<<(std::ostream& os, const GeneRegion<T, D>& x) noexcept {
  os << x.gene.get()
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
}

// derived comparisons (ContigRegion, ContigRegion)
template <class T1, class T2, class D1, class D2>
inline bool operator>(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return y < x;
}
template <class T1, class T2, class D1, class D2>
inline bool operator<=(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return !(y < x);
}
template <class T1, class T2, class D1, class D2>
inline bool operator>=(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return !(x < y);
}

// derived comparisons (GeneRegion, GeneRegion)
template <class T1, class T2, class D1, class D2>
inline bool operator>(
    const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return y < x;
}
template <class T1, class T2, class D1, class D2>
inline bool operator>=(
    const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return !(x < y);
}
template <class T1, class T2, class D1, class D2>
inline bool operator<=(
    const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return !(y < x);
}

// derived comparisons (KnownGene, GeneRegion)
template <class T, class D>
inline bool operator>(const GeneRegion<T, D>& x, const KnownGene& y) noexcept {
  return y < x;
}
template <class T, class D>
inline bool operator>(const KnownGene& x, const GeneRegion<T, D>& y) noexcept {
  return y < x;
}
template <class T, class D>
inline bool operator<=(const GeneRegion<T, D>& x, const KnownGene& y) noexcept {
  return !(y < x);
}
template <class T, class D>
inline bool operator<=(const KnownGene& x, const GeneRegion<T, D>& y) noexcept {
  return !(y < x);
}
template <class T, class D>
inline bool operator>=(const GeneRegion<T, D>& x, const KnownGene& y) noexcept {
  return !(x < y);
}
template <class T, class D>
inline bool operator>=(const KnownGene& x, const GeneRegion<T, D>& y) noexcept {
  return !(x < y);
}

// derived comparisons (ContigRegion, GeneRegion)
template <class T1, class T2, class D1, class D2>
inline bool operator>(
    const GeneRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return y < x;
}
template <class T1, class T2, class D1, class D2>
inline bool operator>(
    const ContigRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return y < x;
}
template <class T1, class T2, class D1, class D2>
inline bool operator<=(
    const GeneRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return !(y < x);
}
template <class T1, class T2, class D1, class D2>
inline bool operator<=(
    const ContigRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return !(y < x);
}
template <class T1, class T2, class D1, class D2>
inline bool operator>=(
    const GeneRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return !(x < y);
}
template <class T1, class T2, class D1, class D2>
inline bool operator>=(
    const ContigRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return !(x < y);
}


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
        return !(a < b && a.contig.known_contigs == b.contig.known_contigs);
        });
  }

 public:
  ContigRegions(vecT&& x, NoCheckValid)
      : elements_{x},
        contigs_{x.empty() ? nullptr : x[0].contig.known_contigs} { }
  ContigRegions(const vecT& x, NoCheckValid)
      : elements_{x},
        contigs_{x.empty() ? nullptr : x[0].contig.known_contigs} { }
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
  friend inline bool operator==(const ContigRegions& x, const ContigRegions& y) {
    return std::tie(x.contigs_, x.elements_) == std::tie(y.contigs_, y.elements_);
  }
};
}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_REGIONS_HPP
