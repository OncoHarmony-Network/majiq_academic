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

  // access
  const KnownContig& parent() const noexcept { return contig; }

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

// ignore data when determining 'equality'
template <typename T1, typename T2, typename D1, typename D2>
inline bool operator==(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    == std::tie(y.contig, y.coordinates, y.strand);
}
// will have types (i.e. GeneRegion) that expose contig() and strand() instead
template <typename T, typename D, typename U,
         std::enable_if_t<
         detail::has_contig_function<U>::value
         && detail::has_strand_function<U>::value, bool> = true>
inline bool operator==(const ContigRegion<T, D>& x, const U& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    == std::tie(y.contig(), y.coordinates, y.strand());
}
template <typename T, typename D, typename U,
         std::enable_if_t<
         detail::has_contig_function<U>::value
         && detail::has_strand_function<U>::value, bool> = true>
inline bool operator==(const U& x, const ContigRegion<T, D>& y) noexcept {
  return std::tie(x.contig(), x.coordinates, x.strand())
    == std::tie(y.contig, y.coordinates, y.strand);
}

// order regions by genomic position and strand, ignoring data
template <typename T1, typename T2, typename D1, typename D2>
inline bool operator<(
    const ContigRegion<T1, D1>& x, const ContigRegion<T2, D2>& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(y.contig, y.coordinates, y.strand);
}
// will have types (i.e. GeneRegion) that expose contig() and strand() instead
template <typename T, typename D, typename U,
         std::enable_if_t<
         detail::has_contig_function<U>::value
         && detail::has_strand_function<U>::value, bool> = true>
inline bool operator<(const ContigRegion<T, D>& x, const U& y) noexcept {
  return std::tie(x.contig, x.coordinates, x.strand)
    < std::tie(y.contig(), y.coordinates, y.strand());
}
template <typename T, typename D, typename U,
         std::enable_if_t<
         detail::has_contig_function<U>::value
         && detail::has_strand_function<U>::value, bool> = true>
inline bool operator<(const U& x, const ContigRegion<T, D>& y) noexcept {
  return std::tie(x.contig(), x.coordinates, x.strand())
    < std::tie(y.contig, y.coordinates, y.strand);
}

template <typename T, typename U = T,
         std::enable_if_t<
           detail::has_contig_function<T>::value
           && detail::has_contig_function<U>::value
           && detail::has_strand_function<T>::value
           && detail::has_strand_function<U>::value, bool
         > = true>
struct CompareContigStranded {
  inline bool operator()(const T& x, const U& y) noexcept {
    return std::tie(x.contig(), x.coordinates, x.strand())
      < std::tie(y.contig(), y.coordinates, y.strand());
  }
  template <typename IT, typename DT>
  inline bool operator()(const T& x, const ContigRegion<IT, DT>& y) {
    return x < y;
  }
  template <typename IT, typename DT>
  inline bool operator()(const ContigRegion<IT, DT>& x, const T& y) {
    return x < y;
  }
};

// unstranded comparison with contig-coordinates
template <typename T, typename U = T>
struct CompareContigUnstranded {
  inline bool operator()(const T& x, const U& y) noexcept {
    constexpr bool T_field = detail::has_contig_field<T>::value;
    constexpr bool T_function = detail::has_contig_function<T>::value;
    constexpr bool U_field = detail::has_contig_field<U>::value;
    constexpr bool U_function = detail::has_contig_function<U>::value;
    static_assert((T_field || T_function),
        "T does not have contig for unstranded contigregion comparison");
    static_assert((U_field || U_function),
        "U does not have contig for unstranded contigregion comparison");
    if constexpr(T_field && U_field) {
      return std::tie(x.contig, x.coordinates)
        < std::tie(y.contig, y.coordinates);
    } else if constexpr(T_field && U_function) {
      return std::tie(x.contig(), x.coordinates)
        < std::tie(y.contig, y.coordinates);
    } else if constexpr(T_function && U_field) {
      return std::tie(x.contig, x.coordinates)
        < std::tie(y.contig(), y.coordinates);
    } else {
      return std::tie(x.contig(), x.coordinates)
        < std::tie(y.contig(), y.coordinates);
    }
  }
};


}  // namespace detail
}  // namespace majiq


#endif  // MAJIQ_CONTIGREGION_HPP
