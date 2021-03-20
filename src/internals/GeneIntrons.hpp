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

#include "GeneRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "GeneConnections.hpp"
#include "Exons.hpp"


namespace majiq {
struct GeneIntron : public detail::GeneConnection<ClosedInterval> {
 public:
  using BaseT = detail::GeneConnection<ClosedInterval>;

  // constructors
  GeneIntron(
      KnownGene gene, ClosedInterval coordinates, detail::ConnectionData data)
      : BaseT{gene, coordinates, data} { }
  GeneIntron(KnownGene gene, ClosedInterval coordinates,
      bool denovo, bool passed_build, bool simplified)
      : GeneIntron{gene, coordinates,
          detail::ConnectionData{denovo, passed_build, simplified}} { }
  GeneIntron(KnownGene gene, ClosedInterval coordinates)
      : GeneIntron{gene, coordinates, false, false, false} { }
  GeneIntron() : GeneIntron{KnownGene{}, ClosedInterval{}} { }
  GeneIntron(const GeneIntron& x) = default;
  GeneIntron(GeneIntron&& x) = default;
  GeneIntron& operator=(const GeneIntron& x) = default;
  GeneIntron& operator=(GeneIntron&& x) = default;
};


class GeneIntrons : public detail::GeneConnections<GeneIntron, false> {
  using BaseT = detail::GeneConnections<GeneIntron, false>;

 public:
  // NOTE: assumes that connections to connected_exons already defined in x
  GeneIntrons(
      const std::shared_ptr<Genes>& genes, std::vector<GeneIntron>&& x,
      const std::shared_ptr<Exons>& connected_exons)
      : BaseT{genes, std::move(x), connected_exons} { }
  GeneIntrons(
      const std::shared_ptr<Genes>& genes, std::vector<GeneIntron>&& x)
      : GeneIntrons{genes, std::move(x), nullptr} { }

  void connect_exons(const std::shared_ptr<Exons>& exons_ptr) override {
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
   * update connection data using overlapping introns from another source
   */
  void UpdateFlagsFrom(const GeneIntrons& donor) const {
    if (parents() != donor.parents()) {
      throw std::invalid_argument(
          "UpdateFlagsFrom requires introns objects to share genes object");
    }
    // otherwise, do per gene
    for (const auto& gene : *parents()) {
      auto it = begin_parent(gene);
      const auto it_end = end_parent(gene);
      for (auto donor_it = donor.begin_parent(gene);
          donor_it != donor.end_parent(gene); ++donor_it) {
        // get first intron that could overlap
        it = std::find_if(it, it_end,
            [&donor_coords = donor_it->coordinates](const GeneIntron& x) {
            return !IntervalPrecedes(x.coordinates, donor_coords); });
        if (it == it_end) { break; }
        for (auto update_it = it;
            (update_it != it_end
             && !IntervalPrecedes(
               donor_it->coordinates, update_it->coordinates));
            ++update_it) {
          // update_it overlaps with donor
          if (!donor_it->denovo()) {
            update_it->denovo() = false;
          }
          if (!donor_it->simplified()) {
            update_it->simplified() = false;
          }
          if (donor_it->passed_build()) {
            update_it->passed_build() = true;
          }
        }  // loop over our introns that overlap with current donor
      }  // loop over donors
    }  // loop over genes
    return;
  }
  /**
   * get all potential introns from given exons
   *
   * Note that automatically connects introns to the exons
   */
  static GeneIntrons PotentialIntrons(
      const std::shared_ptr<Exons>& exons_ptr,
      bool make_simplified) {
    if (exons_ptr == nullptr) {
      throw std::runtime_error("PotentialIntrons requires non-null exons");
    }
    const Exons& exons = *exons_ptr;
    if (exons.parents() == nullptr) {
      throw std::runtime_error("Exons must have non-null genes as parents");
    }

    std::vector<GeneIntron> result_vec;

    for (const auto& gene : *(exons.parents())) {
      const auto exon_it_end = exons.end_parent(gene);
      auto prev_exon_it = exon_it_end;  // haven't gotten to first full exon
      for (auto exon_it = exons.begin_parent(gene);
          exon_it != exon_it_end; ++exon_it) {
        if (exon_it->is_full_exon()) {
          if (prev_exon_it != exon_it_end) {
            size_t start_exon_idx = prev_exon_it - exons.begin();
            size_t end_exon_idx = exon_it - exons.begin();
            result_vec.emplace_back(
                gene,
                ClosedInterval{
                  prev_exon_it->coordinates.end + 1,
                  exon_it->coordinates.start - 1},
                detail::ConnectionData{
                  true,  // consider denovo for now
                  false,  // hasn't passed build
                  make_simplified,
                  start_exon_idx,
                  end_exon_idx});
          }  // handle intron between prev_exon_it and exon_it
          prev_exon_it = exon_it;  // update prev_exon, wait for next full exon
        }  // conditional -- if new full exon
      }  // loop over exons
    }  // loop over genes

    return GeneIntrons{exons.parents(), std::move(result_vec), exons_ptr};
  }
};
}  // namespace majiq

#endif  // MAJIQ_GENEINTRONS_HPP
