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

 public:
  GeneIntrons(
      const typename BaseT::ParentsPtrT& parents, typename BaseT::vecT&& x)
      : BaseT{parents, std::move(x)} { }

  GeneIntrons FilterPassed(bool keep_annotated, bool discard_denovo) const {
    std::vector<GeneIntron> result_vec;
    std::copy_if(begin(), end(), std::back_inserter(result_vec),
        [keep_annotated, discard_denovo](const GeneIntron& x) {
        return (
            x.denovo()
            ? !discard_denovo && x.passed_build()
            : keep_annotated || x.passed_build());
        });
    return GeneIntrons{parents(), std::move(result_vec)};
  }
  // get all potential introns (i.e. denovo as well) compared to exons
  GeneIntrons PotentialIntrons(const Exons& exons) const {
    std::vector<GeneIntron> result_vec;

    if (parents() != exons.parents()) {
      throw std::invalid_argument(
          "PotentialIntrons exons/introns do not share same Genes");
    } else if (parents() != nullptr) {
      constexpr position_t EMPTY = -2;
      for (const auto& gene : *parents()) {
        auto exon_it = exons.begin_parent(gene);
        const auto exon_it_end = exons.end_parent(gene);
        auto intron_it = begin_parent(gene);
        const auto intron_it_end = end_parent(gene);
        position_t last_end = EMPTY;  // end from last exon
        for (; exon_it != exon_it_end; ++exon_it) {
          if (exon_it->is_full_exon()) {
            if (last_end != EMPTY) {
              // coordinates for new intron
              ClosedInterval coordinates = ClosedInterval{
                  last_end + 1, exon_it->coordinates.start - 1};
              // get intron that overlaps, if any
              intron_it = std::find_if(
                  intron_it, intron_it_end,
                  [&coordinates](const GeneIntron& x) {
                  return !IntervalPrecedes(x.coordinates, coordinates);
                  });
              result_vec.emplace_back(
                  gene,
                  ClosedInterval{last_end + 1, exon_it->coordinates.start - 1},
                  // if overlapping intron from self, copy data otherwise denovo
                  (intron_it != intron_it_end
                   && IntervalIntersects(coordinates, intron_it->coordinates))
                  ? intron_it->data : detail::Connection{true, false, false});
            }
            last_end = exon_it->coordinates.end;
          }
        }  // loop over exons for a gene
      }  // loop over genes
    }  // when the containers aren't empty

    return GeneIntrons{parents(), std::move(result_vec)};
  }
};

}  // namespace majiq

#endif  // MAJIQ_GENEINTRONS_HPP
