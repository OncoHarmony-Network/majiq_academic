/**
 * GeneJunctions.hpp
 *
 * Junctions for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_GENEJUNCTIONS_HPP
#define MAJIQ_GENEJUNCTIONS_HPP

#include <tuple>
#include <vector>
#include <memory>
#include <functional>
#include <utility>
#include <sstream>
#include <boost/functional/hash.hpp>

#include "GeneRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "Connection.hpp"
#include "Exons.hpp"


namespace majiq {

struct GeneJunction
    : public detail::GeneRegion<OpenInterval, detail::Connection> {
 public:
  using BaseT = detail::GeneRegion<OpenInterval, detail::Connection>;
  // access data nicely
  const bool& denovo() const noexcept { return data.denovo; }
  bool& denovo() noexcept { return data.denovo; }
  bool& passed_build() const noexcept { return data.passed_build; }
  bool& simplified() const noexcept { return data.simplified; }
  size_t& start_exon_idx() const noexcept { return data.start_exon_idx; }
  size_t& end_exon_idx() const noexcept { return data.end_exon_idx; }
  const size_t& src_exon_idx() const noexcept {
    return gene.strand() == GeneStrandness::FORWARD
      ? start_exon_idx() : end_exon_idx();
  }
  const size_t& dst_exon_idx() const noexcept {
    return gene.strand() == GeneStrandness::FORWARD
      ? end_exon_idx() : start_exon_idx();
  }
  const bool is_exitron() const noexcept { return data.is_exitronic(); }

  // constructors
  GeneJunction(KnownGene gene, OpenInterval coordinates,
      bool denovo, bool passed_build, bool simplified)
      : BaseT{gene, coordinates,
        detail::Connection{denovo, passed_build, simplified}} {
  }
  GeneJunction(KnownGene gene, OpenInterval coordinates)
      : GeneJunction{gene, coordinates, false, false, false} {
  }
  GeneJunction() : GeneJunction{KnownGene{}, OpenInterval{}} { }
  GeneJunction(const GeneJunction& x) = default;
  GeneJunction(GeneJunction&& x) = default;
  GeneJunction& operator=(const GeneJunction& x) = default;
  GeneJunction& operator=(GeneJunction&& x) = default;
};

// override boost::hash
inline std::size_t hash_value(const GeneJunction& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// override std::hash for GeneJunction
namespace std {
template <> struct hash<majiq::GeneJunction> {
  std::size_t operator()(const majiq::GeneJunction& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
class GeneJunctions : public detail::Regions<GeneJunction, true> {
  using BaseT = detail::Regions<GeneJunction, true>;

 private:
  std::shared_ptr<Exons> connected_exons_;

  /**
   * vector of indexes to exons that match junction starts.
   * Throws exception if unable to find all matches
   */
  std::vector<size_t> start_exon_idx(const Exons& exons) const {
    // initialize temporary container for result
    // (don't do the update in place in case there are disconnected junctions)
    std::vector<size_t> result(size());
    size_t eidx = 0;
    for (size_t jidx = 0; jidx < size(); ++jidx) {
      const GeneJunction& j = (*this)[jidx];
      // skip exons that cannot intersect
      for (;
          eidx < exons.size()
          && (exons[eidx].gene < j.gene
            || (exons[eidx].gene == j.gene
              && IntervalPrecedes(exons[eidx].coordinates, j.coordinates)));
          ++eidx) {
      }  // done skipping exons that precede junction
      if (eidx == exons.size()) {
        std::ostringstream oss;
        oss << "Ran out of exons for junction starts (at idx=" << jidx << ")";
        throw std::runtime_error(oss.str());
      } else if (
          // matches gene
          exons[eidx].gene == j.gene
          // matches coordinates
          && (exons[eidx].coordinates.end == j.coordinates.start
            || exons[eidx].coordinates.contains(j.coordinates.start))) {
        result[jidx] = eidx;
      } else {
        std::ostringstream oss;
        oss << "Unable to match junction start (idx=" << jidx
          << ") to exon (expected at idx=" << eidx << ")";
        throw std::runtime_error(oss.str());
      }
    }
    return result;
  }
  std::vector<size_t> end_exon_idx(const Exons& exons) const {
    std::vector<size_t> result(size());
    for (size_t jidx = 0; jidx < size(); ++jidx) {
      const GeneJunction& j = (*this)[jidx];
      const size_t eidx
        = exons.overlap_lower_bound(j.gene, j.coordinates.end) - exons.begin();
      if (eidx < exons.size() && exons[eidx].gene == j.gene
          && (exons[eidx].coordinates.start == j.coordinates.end
            || exons[eidx].coordinates.contains(j.coordinates.end))) {
        result[jidx] = eidx;
      } else {
        std::ostringstream oss;
        oss << "Unable to match junction end (idx=" << jidx
          << ") to exon (expected at idx=" << eidx << ")";
        throw std::runtime_error(oss.str());
      }
    }
    return result;
  }

 public:
  void connect_exons(const std::shared_ptr<Exons>& exons_ptr) {
    if (exons_ptr == nullptr || exons_ptr == connected_exons_) {
      return;  // don't do anything
    }
    const Exons& exons = *exons_ptr;
    if (parents() != exons.parents()) {
      throw std::invalid_argument(
          "junction/exon genes do not match in connect_exons()");
    }
    const auto start_exons = start_exon_idx(exons);
    const auto end_exons = end_exon_idx(exons);
    // no exception --> valid to connect them now in place
    for (size_t i = 0; i < size(); ++i) {
      (*this)[i].start_exon_idx() = start_exons[i];
      (*this)[i].end_exon_idx() = end_exons[i];
    }
    connected_exons_ = exons_ptr;  // update pointer to connected exons
    return;
  }

  bool is_connected() const { return connected_exons_ != nullptr; }
  const std::shared_ptr<Exons>& connected_exons() const {
    return connected_exons_;
  }

  GeneJunctions(
      const std::shared_ptr<Genes>& genes, std::vector<GeneJunction>&& x,
      const std::shared_ptr<Exons>& connected_exons)
      : BaseT{genes, std::move(x)}, connected_exons_{connected_exons} {
    // NOTE: assumes that connections to connected_exons already defined in x
    if (parents() == nullptr) {
      throw std::invalid_argument("GeneJunctions cannot have null genes");
    }
  }
  GeneJunctions(
      const std::shared_ptr<Genes>& genes, std::vector<GeneJunction>&& x)
      : GeneJunctions{genes, std::move(x), nullptr} { }
};
}  // namespace majiq

#endif  // MAJIQ_GENEJUNCTIONS_HPP
