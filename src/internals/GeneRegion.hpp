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

  // access
  const KnownGene& parent() const noexcept { return gene; }

  // expose contig/strand
  const KnownContig& contig() const { return gene.contig(); }
  const GeneStrandness& strand() const { return gene.strand(); }

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

// equality if gene/coordinates match, ignoring data
template <typename T1, typename T2, typename D1, typename D2>
inline bool operator==(
    const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return std::tie(x.gene, x.coordinates) == std::tie(y.gene, y.coordinates);
}
template <typename T1, typename T2, typename D1, typename D2>
inline bool operator!=(
    const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return !(x == y);
}

// order regions by gene, then coordinates
template <typename T1, typename T2, typename D1, typename D2>
inline bool operator<(
    const GeneRegion<T1, D1>& x, const GeneRegion<T2, D2>& y) noexcept {
  return std::tie(x.gene, x.coordinates) < std::tie(y.gene, y.coordinates);
}
// in case we have types that expose gene() instead of gene
template <typename T, typename D, typename U,
         std::enable_if_t<detail::has_gene_function<U>::value> = true>
inline bool operator<(const GeneRegion<T, D>& x, const U& y) noexcept {
  return std::tie(x.gene, x.coordinates)
    < std::tie(y.gene(), y.coordinates);
}
template <typename U, typename T, typename D,
         std::enable_if_t<detail::has_gene_function<U>::value> = true>
inline bool operator<(const U& x, const GeneRegion<T, D>& y) noexcept {
  return std::tie(x.gene(), x.coordinates)
    < std::tie(y.gene, y.coordinates);
}

}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_GENEREGION_HPP
