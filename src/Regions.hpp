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
  ContigRegion(const ContigRegion& x) = default;
  ContigRegion(ContigRegion&& x) = default;
  ContigRegion& operator=(const ContigRegion& x) = default;
  ContigRegion& operator=(ContigRegion&& x) = default;
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
  GeneRegion(const GeneRegion& x) = default;
  GeneRegion(GeneRegion&& x) = default;
  GeneRegion& operator=(const GeneRegion& x) = default;
  GeneRegion& operator=(GeneRegion&& x) = default;
};

// order regions by genomic position and strand
template <class T1, class T2>
inline bool operator<(
    const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
template <class T1, class T2>
inline bool operator==(
    const ContigRegion<T1>& x, const ContigRegion<T2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    == std::tie(y.contig, y.coordinates, y.strand);
}
/**
 * Order regions by gene, then position
 */
template <class T1, class T2>
inline bool operator<(
    const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return std::tie(x.gene, x.coordinates) < std::tie(y.gene, y.coordinates);
}
template <class T1, class T2>
inline bool operator==(
    const GeneRegion<T1>& x, const GeneRegion<T2>& y) noexcept {
  return std::tie(x.gene, x.coordinates) == std::tie(y.gene, y.coordinates);
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


template <
  class RegionT,
  class CompareRegionT = std::less<RegionT>>
class GeneRegions {
 public:
  using IntervalT = typename RegionT::IntervalT;
  using BaseRegion = GeneRegion<IntervalT>;
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

  static bool is_valid(const vecT& x, const std::shared_ptr<Genes>& genes) {
    // 0 elements
    if (x.empty()) {
      return true;
    }
    // first element
    if (genes != x[0].gene.known_genes) {
      return false;
    }
    // subsequent elements
    for (size_t i = 1; i < x.size(); ++i) {
      // if not sorted or if there are equal regions, it is invalid
      if (CompareRegionT{}(x[i], x[i - 1])
          || (static_cast<BaseRegion>(x[i])
            == static_cast<BaseRegion>(x[i - 1]))
          || genes != x[i].gene.known_genes) {
        return false;
      }
    }
    // passed all conditions for validity
    return true;
  }

 public:
  GeneRegions(vecT&& x, NoCheckValid)
      : elements_{x}, genes_{x.empty() ? nullptr : x[0].gene.known_genes} { }
  GeneRegions(const vecT& x, NoCheckValid)
      : elements_{x}, genes_{x.empty() ? nullptr : x[0].gene.known_genes} { }
  explicit GeneRegions(vecT&& x) : GeneRegions{x, NoCheckValid{}} {
    if (!is_valid(x, genes_)) {
      throw std::invalid_argument("vector input to GeneRegions is invalid");
    }
  }
  explicit GeneRegions(const vecT& x) : GeneRegions{x, NoCheckValid{}} {
    if (!is_valid(x, genes_)) {
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
}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_REGIONS_HPP
