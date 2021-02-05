/**
 * ContigRegion.hpp
 *
 * Base class for region sitting on Contig at some interval on one strand vs
 * another
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_CONTIGREGION_HPP
#define MAJIQ_CONTIGREGION_HPP

#include <tuple>

#include "Interval.hpp"
#include "Contigs.hpp"
#include "MajiqTypes.hpp"


namespace majiq {
namespace detail {

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
      KnownContig _contig, IntervalT _coordinates, GeneStrandness _strand, DataT _data)
      : contig{_contig},
        coordinates{_coordinates},
        strand{_strand},
        data{_data} { }
  ContigRegion(KnownContig _contig, IntervalT _coordinates, GeneStrandness _strand)
      : ContigRegion{_contig, _coordinates, _strand, DataT{}} { }
  ContigRegion()
      : ContigRegion{KnownContig{}, IntervalT{}, GeneStrandness::AMBIGUOUS, DataT{}} { }
  ContigRegion(const ContigRegion&) = default;
  ContigRegion(ContigRegion&&) = default;
  ContigRegion& operator=(const ContigRegion&) = default;
  ContigRegion& operator=(ContigRegion&&) = default;

  // drop data between templates
  template <class OtherDataT>
  ContigRegion(const ContigRegion<IntervalT, OtherDataT>& x)
      : ContigRegion{x.contig, x.coordinates, x.strand} { }
  template <class OtherDataT>
  ContigRegion(ContigRegion<IntervalT, OtherDataT>&& x)
      : ContigRegion{x.contig, x.coordinates, x.strand} { }
};

// order regions by genomic position and strand, ignoring data.
template <class T1, class T2, class D1, class D2>
inline bool operator<(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
// ignore data when determining 'equality'
template <class T1, class T2, class D1, class D2>
inline bool operator==(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    == std::tie(y.contig, y.coordinates, y.strand);
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

}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_CONTIGREGION_HPP
