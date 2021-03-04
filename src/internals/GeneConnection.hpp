/**
 * GeneConnection.hpp
 *
 * GeneRegion specialization using ConnectionData
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_GENECONNECTION_HPP
#define MAJIQ_GENECONNECTION_HPP

#include <tuple>

#include "ConnectionData.hpp"
#include "GeneRegion.hpp"
#include "MajiqTypes.hpp"


namespace majiq {
namespace detail {

template <typename IntervalT>
struct GeneConnection
    : public detail::GeneRegion<IntervalT, detail::ConnectionData> {
 public:
  using BaseT = detail::GeneRegion<IntervalT, detail::ConnectionData>;
  using DataT = detail::ConnectionData;
  // constructors
  GeneConnection(KnownGene _gene, IntervalT _coordinates, DataT _data)
      : BaseT{_gene, _coordinates, _data} { }
  GeneConnection(KnownGene _gene, IntervalT _coordinates)
      : GeneConnection{_gene, _coordinates, DataT{}} {
  }
  GeneConnection() : GeneConnection{KnownGene{}, IntervalT{}} { }
  GeneConnection(const GeneConnection&) = default;
  GeneConnection(GeneConnection&&) = default;
  GeneConnection& operator=(const GeneConnection&) = default;
  GeneConnection& operator=(GeneConnection&&) = default;

  // access data nicely
  const bool& denovo() const noexcept { return this->data.denovo; }
  bool& denovo() noexcept { return this->data.denovo; }
  bool& passed_build() const noexcept { return this->data.passed_build; }
  bool& simplified() const noexcept { return this->data.simplified; }
  size_t& start_exon_idx() const noexcept { return this->data.start_exon_idx; }
  size_t& end_exon_idx() const noexcept { return this->data.end_exon_idx; }
  const size_t& src_exon_idx() const noexcept {
    return this->gene.strand() == GeneStrandness::FORWARD
      ? start_exon_idx() : end_exon_idx();
  }
  const size_t& dst_exon_idx() const noexcept {
    return this->gene.strand() == GeneStrandness::FORWARD
      ? end_exon_idx() : start_exon_idx();
  }
  const size_t& ref_exon_idx(const EventType& type) const noexcept {
    return type == EventType::SRC_EVENT ? src_exon_idx() : dst_exon_idx();
  }
  const size_t& other_exon_idx(const EventType& type) const noexcept {
    return type == EventType::SRC_EVENT ? dst_exon_idx() : src_exon_idx();
  }
  // first element is reference coordinate, second is other coordinate
  std::tuple<const position_t&, const position_t&>
  coordinates_ref_other(const EventType& type) const noexcept {
    return this->strand_forward() == (type == EventType::SRC_EVENT)
      ? this->coordinates.as_tuple() : this->coordinates.rev_tuple();
  }
  const position_t& ref_coordinate(const EventType& type) const noexcept {
    return std::get<0>(coordinates_ref_other(type));
  }
  const position_t& other_coordinate(const EventType& type) const noexcept {
    return std::get<1>(coordinates_ref_other(type));
  }

  // helpers
  bool is_exitron() const noexcept { return this->data.is_exitronic(); }
  bool for_event() const noexcept { return !(is_exitron() || simplified()); }
  bool for_passed() const noexcept { return passed_build() && for_event(); }
};

}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_GENECONNECTION_HPP
