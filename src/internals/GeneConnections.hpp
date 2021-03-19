/**
 * GeneConnections.hpp
 *
 * GeneRegion specialization using ConnectionData, container over them
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_GENECONNECTIONS_HPP
#define MAJIQ_GENECONNECTIONS_HPP

#include <tuple>

#include "ConnectionData.hpp"
#include "GeneRegion.hpp"
#include "MajiqTypes.hpp"
#include "Exons.hpp"
#include "checksum.hpp"


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
  bool& denovo() const noexcept { return this->data.denovo; }
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

template <typename GeneConnectionT, bool HAS_OVERLAPS>
class GeneConnections : public Regions<GeneConnectionT, HAS_OVERLAPS> {
  using BaseT = detail::Regions<GeneConnectionT, HAS_OVERLAPS>;

 protected:
  std::shared_ptr<Exons> connected_exons_;

 public:
  GeneConnections(
      const std::shared_ptr<Genes>& genes, std::vector<GeneConnectionT>&& x,
      const std::shared_ptr<Exons>& connected_exons)
      : BaseT{genes, std::move(x)}, connected_exons_{connected_exons} {
    if (this->parents() == nullptr) {
      throw std::invalid_argument("GeneConnections cannot have null genes");
    }
  }

  bool is_connected() const { return connected_exons_ != nullptr; }
  const std::shared_ptr<Exons>& connected_exons() const {
    return connected_exons_;
  }
  // connect exons to held GeneConnectionT values. Sets connected_exons_
  virtual void connect_exons(const std::shared_ptr<Exons>& exons_ptr) = 0;

  void pass_all() const {
    std::for_each(this->begin(), this->end(),
        [](const GeneConnectionT& x) { x.passed_build() = true; });
  }
  void simplify_all() const {
    std::for_each(this->begin(), this->end(),
        [](const GeneConnectionT& x) { x.simplified() = true; });
  }
  void unsimplify_all() const {
    std::for_each(this->begin(), this->end(),
        [](const GeneConnectionT& x) { x.simplified() = false; });
  }

  template <bool PROCESS_DATA>
  checksum_t checksum() const {
    checksum_gen_t gen;
    for (const auto& x : *this) {
      gen.process_bytes(&x.gene.idx_, sizeof(x.gene.idx_));
      gen.process_bytes(&x.coordinates.start, sizeof(x.coordinates.start));
      gen.process_bytes(&x.coordinates.end, sizeof(x.coordinates.end));
      if constexpr(PROCESS_DATA) {
        gen.process_bytes(&x.data.denovo, sizeof(x.data.denovo));
        gen.process_bytes(&x.data.passed_build, sizeof(x.data.passed_build));
        gen.process_bytes(&x.data.simplified, sizeof(x.data.simplified));
        gen.process_bytes(
            &x.data.start_exon_idx, sizeof(x.data.start_exon_idx));
        gen.process_bytes(&x.data.end_exon_idx, sizeof(x.data.end_exon_idx));
      }
    }
    return checksum_t{gen.checksum()};
  }
};

}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_GENECONNECTIONS_HPP
