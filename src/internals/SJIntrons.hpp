/**
 * SJIntrons.hpp
 *
 * Introns on contigs for quantification
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_SJINTRONS_HPP
#define MAJIQ_SJINTRONS_HPP

#include <vector>
#include <memory>
#include <utility>

#include "ContigRegion.hpp"
#include "Regions.hpp"
#include "Exons.hpp"
#include "GeneIntrons.hpp"

namespace majiq {

namespace detail {
struct AnnotatedIntronStatus {
  bool annotated_;
  explicit AnnotatedIntronStatus(bool annotated) : annotated_{annotated} { }
  AnnotatedIntronStatus() : AnnotatedIntronStatus{false} { }
  AnnotatedIntronStatus(const AnnotatedIntronStatus&) = default;
  AnnotatedIntronStatus(AnnotatedIntronStatus&&) = default;
  AnnotatedIntronStatus& operator=(const AnnotatedIntronStatus&) = default;
  AnnotatedIntronStatus& operator=(AnnotatedIntronStatus&&) = default;
};
}  // namespace detail

struct SJIntron
    : public detail::ContigRegion<
        ClosedInterval, detail::AnnotatedIntronStatus> {
 public:
  using AnnotatedIntronStatus = detail::AnnotatedIntronStatus;
  using BaseT = detail::ContigRegion<ClosedInterval, AnnotatedIntronStatus>;
  const bool& annotated() const noexcept { return data.annotated_; }
  bool& annotated() noexcept { return data.annotated_; }

  SJIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand, AnnotatedIntronStatus annotated)
      : BaseT{contig, coordinates, strand, annotated} { }
  SJIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand, bool annotated)
      : SJIntron{contig, coordinates, strand,
        AnnotatedIntronStatus{annotated}} { }
  SJIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand)
      : SJIntron{contig, coordinates, strand, AnnotatedIntronStatus{}} { }
  SJIntron()
      : SJIntron{KnownContig{}, ClosedInterval{},
        GeneStrandness::AMBIGUOUS} { }
  SJIntron(const SJIntron&) = default;
  SJIntron(SJIntron&&) = default;
  SJIntron& operator=(const SJIntron&) = default;
  SJIntron& operator=(SJIntron&&) = default;
};

class SJIntrons : public detail::Regions<SJIntron, true> {
  using BaseT = detail::Regions<SJIntron, true>;
 public:
  static SJIntrons FromGeneExonsAndIntrons(
      const Exons& exons, const GeneIntrons& gene_introns, const bool stranded);

  SJIntrons(
      const std::shared_ptr<Contigs>& contigs, std::vector<SJIntron>&& x)
      : BaseT{contigs, std::move(x)} { }
};

}  // namespace majiq

#endif
