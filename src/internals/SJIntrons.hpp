/**
 * SJIntrons.hpp
 *
 * Introns on contigs for quantification
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_SJINTRONS_HPP
#define MAJIQ_SJINTRONS_HPP

#include <algorithm>
#include <vector>
#include <map>
#include <memory>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "ContigRegion.hpp"
#include "Regions.hpp"
#include "Exons.hpp"
#include "GeneIntrons.hpp"

namespace majiq {

// data for SJ introns (i.e. if for annotated gene introns only vs not)
struct SJIntronStatus {
  bool annotated_;
  explicit SJIntronStatus(bool annotated) : annotated_{annotated} { }
  SJIntronStatus() : SJIntronStatus{false} { }
  SJIntronStatus(const SJIntronStatus&) = default;
  SJIntronStatus(SJIntronStatus&&) = default;
  SJIntronStatus& operator=(const SJIntronStatus&) = default;
  SJIntronStatus& operator=(SJIntronStatus&&) = default;
};

struct SJIntron
    : public detail::ContigRegion<ClosedInterval, SJIntronStatus> {
 public:
  using BaseT = detail::ContigRegion<ClosedInterval, SJIntronStatus>;
  const bool& annotated() const noexcept { return data.annotated_; }
  bool& annotated() noexcept { return data.annotated_; }

  SJIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand, SJIntronStatus annotated)
      : BaseT{contig, coordinates, strand, annotated} { }
  SJIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand, bool annotated)
      : SJIntron{contig, coordinates, strand,
        SJIntronStatus{annotated}} { }
  SJIntron(KnownContig contig, ClosedInterval coordinates,
      GeneStrandness strand)
      : SJIntron{contig, coordinates, strand, SJIntronStatus{}} { }
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
  inline static SJIntrons FromGeneExonsAndIntrons(
      const Exons& exons, const GeneIntrons& gene_introns, const bool stranded);

  SJIntrons(
      const std::shared_ptr<Contigs>& contigs, std::vector<SJIntron>&& x)
      : BaseT{contigs, std::move(x)} { }
};

namespace detail {

// types over gene intron/exon coordinates for desired ordering
enum class SJIntronEvidenceType : unsigned char {
  FIRST_EXON_START,
  EXON_START,
  ANNOTATED_INTRON_END,
  ANNOTATED_INTRON_START,
  EXON_END,
  LAST_EXON_END
};

// this evidence is sorted after position
using SJIntronEvidence = std::pair<position_t, SJIntronEvidenceType>;

/**
 * Given overlapping genes from [first, last), accumulate evidence for introns
 */
inline std::map<GeneStrandness, std::vector<SJIntronEvidence>>
SJIntronEvidenceForOverGene(
    const Exons& exons, const GeneIntrons& gene_introns,
    const KnownGene& first, const KnownGene& last, bool stranded) {
  std::map<GeneStrandness, std::vector<SJIntronEvidence>> result;
  for (KnownGene gene = first; gene < last; ++gene) {
    // evidence from gene's exons
    for (auto it = exons.begin_parent(gene);
        it != exons.end_parent(gene); ++it) {
      if (it->coordinates.is_full_interval()) {
        GeneStrandness strand
          = stranded ? it->gene.strand() : GeneStrandness::AMBIGUOUS;
        result[strand].emplace_back(it->coordinates.start,
            it == exons.begin_parent(gene)
            ? SJIntronEvidenceType::FIRST_EXON_START
            : SJIntronEvidenceType::EXON_START);
        result[strand].emplace_back(it->coordinates.end,
            it == std::prev(exons.end_parent(gene))
            ? SJIntronEvidenceType::LAST_EXON_END
            : SJIntronEvidenceType::EXON_END);
      }
    }
    // evidence from gene's introns
    for (auto it = gene_introns.begin_parent(gene);
        it != gene_introns.end_parent(gene); ++it) {
      if (!(it->denovo()) && it->coordinates.is_full_interval()) {
        GeneStrandness strand
          = stranded ? it->gene.strand() : GeneStrandness::AMBIGUOUS;
        // NOTE: we adjust intron coordinates +/- 1 back to exon boundaries
        result[strand].emplace_back(it->coordinates.start - 1,
            SJIntronEvidenceType::ANNOTATED_INTRON_START);
        result[strand].emplace_back(it->coordinates.end + 1,
            SJIntronEvidenceType::ANNOTATED_INTRON_END);
      }
    }
  }
  // sort evidence for each strand
  for (auto&& [strand, evidence_vec] : result) {
    std::sort(evidence_vec.begin(), evidence_vec.end());
  }
  return result;
}

/**
 * Add introns inferred from evidence into result
 */
inline void AddIntronsFromEvidence(KnownContig contig, GeneStrandness strand,
    const std::vector<SJIntronEvidence>& evidence,
    std::vector<SJIntron>& result) {
  int exon_ct = 0;
  int annotated_ct = 0;
  int intron_ct = 0;
  constexpr position_t EMPTY = -2;
  position_t intron_start = EMPTY;
  for (const auto& [position, from] : evidence) {
    switch (from) {
      case SJIntronEvidenceType::EXON_START:
        --intron_ct;
        // no break, do everything in FIRST_EXON_START too
        // (i.e. FIRST_EXON_START just doesn't end an intron)
      case SJIntronEvidenceType::FIRST_EXON_START:
        ++exon_ct;
        if (intron_start != EMPTY) {
          result.emplace_back(
              contig, ClosedInterval{intron_start, position - 1}, strand,
              annotated_ct > 0 /* annotated? */);
          intron_start = EMPTY;
        }
        break;
      case SJIntronEvidenceType::EXON_END:
        ++intron_ct;
        // no break, do everything in LAST_EXON_END too
        // (i.e. LAST_EXON_END just doesn't start an intron)
      case SJIntronEvidenceType::LAST_EXON_END:
        --exon_ct;
        if (intron_ct > 0 && exon_ct == 0) {
          intron_start = position + 1;
        }
        break;
      case SJIntronEvidenceType::ANNOTATED_INTRON_END:
        --annotated_ct;
        break;
      case SJIntronEvidenceType::ANNOTATED_INTRON_START:
        ++annotated_ct;
        break;
    }
  }
}

}  // namespace detail

// implementation for SJIntrons::FromGeneExonsAndIntrons
SJIntrons SJIntrons::FromGeneExonsAndIntrons(
    const Exons& exons, const GeneIntrons& gene_introns, const bool stranded) {
  if (exons.parents_ != gene_introns.parents_) {
    throw std::invalid_argument(
        "SJIntrons gene exons and introns do not share same genes");
  }

  std::vector<SJIntron> result_vec;
  // if there are no exons, there can be no introns
  if (exons.empty()) {
    return SJIntrons{exons.parents()->parents(), std::move(result_vec)};
  }

  // otherwise, operate on sets of genes at a time that overlap
  const auto& genes_ptr = exons.parents_;
  KnownGene gene_it = genes_ptr->begin();
  for (KnownGene gene_it = genes_ptr->begin(),
      next_it = gene_it.NextOverGeneStart();
      gene_it != genes_ptr->end();
      gene_it = next_it, next_it = next_it.NextOverGeneStart()) {
    // get evidence from current set of overlapping genes
    auto stranded_evidence = detail::SJIntronEvidenceForOverGene(
        exons, gene_introns, gene_it, next_it, stranded);
    // update result_vec using evidence
    // each strand gets put in in sorted order, so we can set up merges between
    // them if necessary.
    // NOTE: we only expect to ever see two outputs. Three would be unusual
    // (suggesting third strand type, which is impossible), but we handle this
    // impossibility rather than checking for it and throwing an error
    std::vector<size_t> prev_sizes;
    prev_sizes.reserve(1 + stranded_evidence.size());
    prev_sizes.push_back(result_vec.size());
    for (const auto& [strand, evidence] : stranded_evidence) {
      detail::AddIntronsFromEvidence(
          gene_it.contig(), strand, evidence, result_vec);
      prev_sizes.push_back(result_vec.size());
    }
    for (size_t merge_idx = 2; merge_idx < prev_sizes.size(); ++merge_idx) {
      std::inplace_merge(
          result_vec.begin() + prev_sizes[0],
          result_vec.begin() + prev_sizes[merge_idx - 1],
          result_vec.begin() + prev_sizes[merge_idx]);
    }
  }
  // get final result
  return SJIntrons{exons.parents()->parents(), std::move(result_vec)};
}

}  // namespace majiq

#endif
