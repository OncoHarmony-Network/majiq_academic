/**
 * GeneIntrons.hpp
 *
 * Gene introns for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_GENEINTRONS_HPP
#define MAJIQ_GENEINTRONS_HPP

#include <tuple>
#include <utility>
#include <functional>
#include <algorithm>
#include <vector>
#include <memory>
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
struct GeneIntron
    : public detail::GeneRegion<ClosedInterval, detail::Connection> {
 public:
  using BaseT = detail::GeneRegion<ClosedInterval, detail::Connection>;
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
  GeneIntron(
      KnownGene gene, ClosedInterval coordinates, detail::Connection data)
      : BaseT{gene, coordinates, data} { }
  GeneIntron(KnownGene gene, ClosedInterval coordinates,
      bool denovo, bool passed_build, bool simplified)
      : GeneIntron{gene, coordinates,
          detail::Connection{denovo, passed_build, simplified}} { }
  GeneIntron(KnownGene gene, ClosedInterval coordinates)
      : GeneIntron{gene, coordinates, false, false, false} { }
  GeneIntron() : GeneIntron{KnownGene{}, ClosedInterval{}} { }
  GeneIntron(const GeneIntron& x) = default;
  GeneIntron(GeneIntron&& x) = default;
  GeneIntron& operator=(const GeneIntron& x) = default;
  GeneIntron& operator=(GeneIntron&& x) = default;
};

// override boost hashing
inline std::size_t hash_value(const GeneIntron& x) noexcept {
  std::size_t result = hash_value(x.gene);
  boost::hash_combine(result, x.coordinates);
  return result;
}
}  // namespace majiq
// specialize std::hash for GeneIntron
namespace std {
template <> struct hash<majiq::GeneIntron> {
  std::size_t operator()(const majiq::GeneIntron& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {

class GeneIntrons : public detail::Regions<GeneIntron, false> {
  using BaseT = detail::Regions<GeneIntron, false>;

 private:
  std::shared_ptr<Exons> connected_exons_;

 public:
  GeneIntrons(
      const std::shared_ptr<Genes>& genes, std::vector<GeneIntron>&& x,
      const std::shared_ptr<Exons>& connected_exons)
      : BaseT{genes, std::move(x)}, connected_exons_{connected_exons} {
    // NOTE: assumes that connections to connected_exons already defined in x
    if (parents() == nullptr) {
      throw std::invalid_argument("GeneIntrons cannot have null genes");
    }
  }
  GeneIntrons(
      const std::shared_ptr<Genes>& genes, std::vector<GeneIntron>&& x)
      : GeneIntrons{genes, std::move(x), nullptr} { }

  bool is_connected() const { return connected_exons_ != nullptr; }
  const std::shared_ptr<Exons>& connected_exons() const {
    return connected_exons_;
  }

  void connect_exons(const std::shared_ptr<Exons>& exons_ptr) {
    if (exons_ptr == nullptr || exons_ptr == connected_exons_) {
      return;  // don't do anything
    }
    const Exons& exons = *exons_ptr;
    if (parents() != exons.parents()) {
      throw std::invalid_argument(
          "junction/exon genes do not match in connect_exons()");
    }
    struct _ExonIndexes {
      size_t start;
      size_t end;
    };
    std::vector<_ExonIndexes> result(size());
    for (const auto& gene : *parents()) {
      auto exon_it = exons.begin_parent(gene);
      const auto exon_it_end = exons.end_parent(gene);
      for (auto intron_it = begin_parent(gene);
          intron_it != end_parent(gene);
          ++intron_it) {
        _ExonIndexes intron_exons{};  // get exons for this intron
        exon_it = std::find_if(exon_it, exon_it_end,
            [&intron = *intron_it](const Exon& x) {
            return 1 + x.coordinates.end >= intron.coordinates.start; });
        if (exon_it == exon_it_end
            || !exon_it->is_full_exon()
            || 1 + exon_it->coordinates.end != intron_it->coordinates.start) {
          std::ostringstream oss;
          oss << "Cannot match intron (idx=" << (intron_it - begin())
            << ") to appropriate start exon (expected idx="
            << (exon_it - exons.begin()) << ")";
          throw std::runtime_error(oss.str());
        } else {
          intron_exons.start = exon_it - exons.begin();
        }
        // should be connected to next full exon
        exon_it = std::find_if(1 + exon_it, exon_it_end,
            [](const Exon& x) { return x.is_full_exon(); });
        if (exon_it == exon_it_end
            || 1 + intron_it->coordinates.end != exon_it->coordinates.start) {
          std::ostringstream oss;
          oss << "Cannot match intron (idx=" << (intron_it - begin())
            << ") to appropriate end exon (expected idx="
            << (exon_it - exons.begin()) << ")";
          throw std::runtime_error(oss.str());
        } else {
          intron_exons.end = exon_it - exons.begin();
        }
        result[intron_it - begin()] = intron_exons;
      }  // done looping over introns for gene
    }  // done looping over genes
    // no exceptions -> exons are compatible with introns, update in place now
    for (size_t i = 0; i < size(); ++i) {
      (*this)[i].start_exon_idx() = result[i].start;
      (*this)[i].end_exon_idx() = result[i].end;
    }
    connected_exons_ = exons_ptr;  // update pointer to connected exons
    return;
  }

  GeneIntrons FilterPassed(bool keep_annotated, bool discard_denovo) const {
    std::vector<GeneIntron> result_vec;
    std::copy_if(begin(), end(), std::back_inserter(result_vec),
        [keep_annotated, discard_denovo](const GeneIntron& x) {
        return (
            x.denovo()
            ? !discard_denovo && x.passed_build()
            : keep_annotated || x.passed_build());
        });
    return GeneIntrons{parents(), std::move(result_vec), connected_exons_};
  }
  /**
   * get all potential introns (i.e. denovo as well) compared to exons
   *
   * Note that automatically connects introns to the exons
   */
  GeneIntrons PotentialIntrons(const std::shared_ptr<Exons>& exons_ptr) const {
    if (exons_ptr == nullptr) {
      throw std::runtime_error("PotentialIntrons requires non-null exons");
    }
    const Exons& exons = *exons_ptr;
    std::vector<GeneIntron> result_vec;

    if (parents() != exons.parents()) {
      throw std::invalid_argument(
          "PotentialIntrons exons/introns do not share same Genes");
    } else if (parents() != nullptr) {
      for (const auto& gene : *parents()) {
        auto exon_it = exons.begin_parent(gene);
        const auto exon_it_end = exons.end_parent(gene);
        auto intron_it = begin_parent(gene);
        const auto intron_it_end = end_parent(gene);
        auto prev_exon_it = exon_it_end;  // track previous exon
        for (; exon_it != exon_it_end; ++exon_it) {
          if (exon_it->is_full_exon()) {
            if (prev_exon_it != exon_it_end) {
              // coordinates for new intron
              ClosedInterval coordinates = ClosedInterval{
                  prev_exon_it->coordinates.end + 1,
                  exon_it->coordinates.start - 1};
              // get intron that overlaps, if any
              intron_it = std::find_if(
                  intron_it, intron_it_end,
                  [&coordinates](const GeneIntron& x) {
                  return !IntervalPrecedes(x.coordinates, coordinates);
                  });
              // indexes of exons connected to
              size_t start_idx = prev_exon_it - exons.begin();
              size_t end_idx = exon_it - exons.begin();
              result_vec.emplace_back(
                  gene,
                  ClosedInterval{
                    prev_exon_it->coordinates.end + 1,
                    exon_it->coordinates.start - 1},
                  // if overlapping intron from self, copy data otherwise denovo
                  (intron_it != intron_it_end
                   && IntervalIntersects(coordinates, intron_it->coordinates))
                  ? detail::Connection{intron_it->data, start_idx, end_idx}
                  : detail::Connection{true, false, false, start_idx, end_idx});
            }
            prev_exon_it = exon_it;
          }
        }  // loop over exons for a gene
      }  // loop over genes
    }  // when the containers aren't empty

    return GeneIntrons{parents(), std::move(result_vec), exons_ptr};
  }
};

}  // namespace majiq

#endif  // MAJIQ_GENEINTRONS_HPP
