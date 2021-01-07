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
#include "utils/ContiguousSet.hpp"


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


template <
  class RegionT,
  class CompareRegionT = std::less<RegionT>,
  class SetUnion = utils::DefaultSetUnion<RegionT, CompareRegionT>>
class GeneRegions
    : public utils::contiguous_set<RegionT, CompareRegionT, SetUnion> {
 public:
  using IntervalT = typename RegionT::IntervalT;
  static_assert(std::is_base_of<GeneRegion<IntervalT>, RegionT>::value,
      "GeneRegions type must be subclass of GeneRegion");

  using BaseT = utils::contiguous_set<RegionT, CompareRegionT, SetUnion>;
  using typename BaseT::vecT;
  using BaseT::is_contiguous;
  using BaseT::make_contiguous;
  using BaseT::size;
  using BaseT::sorted_vec_;
  using BaseT::lazy_set_;

  void remap_genes(const std::shared_ptr<Genes>& new_genes) {
    // remap genes
    if (is_contiguous()) {
      // remap in place
      std::for_each(sorted_vec_.begin(), sorted_vec_.end(),
          [&new_genes](RegionT& x) { x.gene = x.gene.remapped(new_genes); });
    } else {
      vecT temp;
      temp.reserve(size());
      auto remapper = [&new_genes](RegionT x) {
        x.gene = x.gene.remapped(new_genes);
        return x;
      };
      std::transform(sorted_vec_.begin(), sorted_vec_.end(),
          std::back_inserter(temp), remapper);
      std::transform(lazy_set_.begin(), lazy_set_.end(),
          std::back_inserter(temp), remapper);
      sorted_vec_ = std::move(temp);
      lazy_set_.clear();
    }
    std::sort(sorted_vec_.begin(), sorted_vec_.end());
  }

  // constructors
  GeneRegions() = default;
  GeneRegions(const GeneRegions& x) = default;
  GeneRegions(GeneRegions&& x) = default;
  GeneRegions& operator=(const GeneRegions& x) = default;
  GeneRegions& operator=(GeneRegions&& x) = default;
  // construct and remap
  GeneRegions(const std::shared_ptr<Genes>& new_genes, GeneRegions&& x)
      : GeneRegions{x} {
    remap_genes(new_genes);
  }
  GeneRegions(
      const std::shared_ptr<Genes>& new_genes,
      std::initializer_list<const GeneRegions&> regions_list)
      : GeneRegions{} {
    for (const auto& regions : regions_list) {
      auto remap_and_insert = [this, &new_genes](RegionT x) {
        x.gene = x.gene.remapped(new_genes);
        insert(x);
        return;
      };
      for_each(regions.sorted_vec_.begin(), regions.sorted_vec_.end(),
          remap_and_insert);
      for_each(regions.lazy_set_.begin(), regions.lazy_set_.end(),
          remap_and_insert);
    }
    make_contiguous();
  }
};
}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_REGIONS_HPP
