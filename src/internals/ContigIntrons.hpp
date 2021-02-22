/**
 * ContigIntrons.hpp
 *
 * Introns on contigs for quantification
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_CONTIGINTRONS_HPP
#define MAJIQ_CONTIGINTRONS_HPP

#include <vector>

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

struct ContigIntron
    : public detail::ContigRegion<
        ClosedInterval, detail::AnnotatedIntronStatus> {
 public:
  using AnnotatedIntronStatus = detail::AnnotatedIntronStatus;
  using BaseT = detail::ContigRegion<ClosedInterval, AnnotatedIntronStatus>;
  const bool& annotated() const noexcept { return data.annotated_; }
  bool& annotated() noexcept { return data.annotated_; }

  ContigIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand, AnnotatedIntronStatus annotated)
      : BaseT{contig, coordinates, strand, annotated} { }
  ContigIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand, bool annotated)
      : ContigIntron{contig, coordinates, strand,
        AnnotatedIntronStatus{annotated}} { }
  ContigIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand)
      : ContigIntron{contig, coordinates, strand, AnnotatedIntronStatus{}} { }
  ContigIntron()
      : ContigIntron{KnownContig{}, ClosedInterval{},
        GeneStrandness::AMBIGUOUS} { }
  ContigIntron(const ContigIntron&) = default;
  ContigIntron(ContigIntron&&) = default;
  ContigIntron& operator=(const ContigIntron&) = default;
  ContigIntron& operator=(ContigIntron&&) = default;
};

class ContigIntrons : public detail::Regions<ContigIntron, true> {
  using BaseT = detail::Regions<ContigIntron, true>;
 public:
  static ContigIntrons FromGeneExonsAndIntrons(
      const Exons& exons, const GeneIntrons& gene_introns, const bool stranded);

  explicit ContigIntrons(std::vector<ContigIntron>&& x)
      : BaseT{std::move(x)} { }
  ContigIntrons() : BaseT{} { }
  ContigIntrons(const ContigIntrons&) = default;
  ContigIntrons(ContigIntrons&&) = default;
  ContigIntrons& operator=(const ContigIntrons&) = default;
  ContigIntrons& operator=(ContigIntrons&&) = default;
};

}  // namespace majiq

#endif
