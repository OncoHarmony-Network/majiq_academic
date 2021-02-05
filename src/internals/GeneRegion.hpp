/**
 * GeneRegion.hpp
 *
 * Base class for intervals belonging to specific genes
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_GENEREGION_HPP
#define MAJIQ_GENEREGION_HPP

#include <tuple>

#include "Interval.hpp"
#include "ContigRegion.hpp"
#include "Genes.hpp"
#include "MajiqTypes.hpp"

namespace majiq {
namespace detail {
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
  GeneRegion(const GeneRegion&) = default;
  GeneRegion(GeneRegion&&) = default;
  GeneRegion& operator=(const GeneRegion&) = default;
  GeneRegion& operator=(GeneRegion&&) = default;
};

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
std::ostream& operator<<(std::ostream& os, const GeneRegion<T, D>& x) noexcept {
  os << x.gene.get()
    << ":"<< x.coordinates.start << "-" << x.coordinates.end;
  return os;
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

}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_GENEREGION_HPP
