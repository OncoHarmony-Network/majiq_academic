/**
 * ContigIntrons.cpp
 *
 * Implementation of ContigIntrons
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "ContigIntrons.hpp"

#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "ContigRegion.hpp"


namespace majiq {

// types over gene intron/exon coordinates for desired ordering
enum class ContigIntronEvidenceType : unsigned char {
  FIRST_EXON_START,
  EXON_START,
  ANNOTATED_INTRON_END,
  ANNOTATED_INTRON_START,
  EXON_END,
  LAST_EXON_END
};

using evidence_t = std::pair<position_t, ContigIntronEvidenceType>;

std::map<GeneStrandness, std::vector<evidence_t>> OverGeneEvidence(
    const Exons& exons, const Introns& introns,
    size_t gene_idx, size_t over_idx, bool stranded) {
  std::map<GeneStrandness, std::vector<evidence_t>> result;
  for (; gene_idx < over_idx; ++gene_idx) {
    // evidence from gene's exons
    for (auto it = exons.begin_parent(gene_idx);
        it != exons.end_parent(gene_idx); ++it) {
      if (it->coordinates.is_full_interval()) {
        GeneStrandness strand
          = stranded ? it->gene.strand() : GeneStrandness::AMBIGUOUS;
        result[strand].emplace_back(it->coordinates.start,
            it == exons.begin_parent(gene_idx)
            ? ContigIntronEvidenceType::FIRST_EXON_START
            : ContigIntronEvidenceType::EXON_START);
        result[strand].emplace_back(it->coordinates.end,
            it == std::prev(exons.end_parent(gene_idx))
            ? ContigIntronEvidenceType::LAST_EXON_END
            : ContigIntronEvidenceType::EXON_END);
      }
    }
    // evidence from gene's introns
    for (auto it = introns.begin_parent(gene_idx);
        it != introns.end_parent(gene_idx); ++it) {
      if (!(it->denovo()) && it->coordinates.is_full_interval()) {
        GeneStrandness strand
          = stranded ? it->gene.strand() : GeneStrandness::AMBIGUOUS;
        // NOTE: we adjust intron coordinates +/- 1 back to exon boundaries
        result[strand].emplace_back(it->coordinates.start - 1,
            ContigIntronEvidenceType::ANNOTATED_INTRON_START);
        result[strand].emplace_back(it->coordinates.end + 1,
            ContigIntronEvidenceType::ANNOTATED_INTRON_END);
      }
    }
  }  // done combining information from all overlapping genes
  // sort evidence for each strand
  for (auto&& [strand, evidence_vec] : result) {
    std::sort(evidence_vec.begin(), evidence_vec.end());
  }
  return result;
}

// add introns inferred from evidence into result
void AddIntronsFromEvidence(KnownContig contig, GeneStrandness strand,
    const std::vector<evidence_t>& evidence,
    std::vector<ContigIntron>& result) {
  int exon_ct = 0;
  int annotated_ct = 0;
  int intron_ct = 0;
  constexpr position_t EMPTY = -2;
  position_t intron_start = EMPTY;
  for (const auto& [position, from] : evidence) {
    switch (from) {
      case ContigIntronEvidenceType::EXON_START:
        --intron_ct;
        // no break, do everything in FIRST_EXON_START too
        // (i.e. FIRST_EXON_START just doesn't end an intron)
      case ContigIntronEvidenceType::FIRST_EXON_START:
        ++exon_ct;
        if (intron_start != EMPTY) {
          result.emplace_back(
              contig, ClosedInterval{intron_start, position - 1}, strand,
              annotated_ct > 0 /* annotated? */);
          intron_start = EMPTY;
        }
        break;
      case ContigIntronEvidenceType::EXON_END:
        ++intron_ct;
        // no break, do everything in LAST_EXON_END too
        // (i.e. LAST_EXON_END just doesn't start an intron)
      case ContigIntronEvidenceType::LAST_EXON_END:
        --exon_ct;
        if (intron_ct > 0 && exon_ct == 0) {
          intron_start = position + 1;
        }
        break;
      case ContigIntronEvidenceType::ANNOTATED_INTRON_END:
        --annotated_ct;
        break;
      case ContigIntronEvidenceType::ANNOTATED_INTRON_START:
        ++annotated_ct;
        break;
    }
  }
}

ContigIntrons ContigIntrons::FromGeneExonsAndIntrons(
    const Exons& exons, const Introns& introns, const bool stranded) {
  if (exons.parents_ != introns.parents_) {
    throw std::invalid_argument(
        "ContigIntrons gene exons and introns do not share same genes");
  }
  // if there are no exons, there can be no introns
  if (exons.empty()) { return ContigIntrons{}; }

  std::vector<ContigIntron> result_vec;
  // otherwise, operate on sets of genes at a time that overlap
  const auto& genes_ptr = exons.parents_;
  size_t gene_idx = 0;
  while (gene_idx < genes_ptr->size()) {
    // get first over_idx past end, non-matching contig, or start past end
    size_t over_idx = 1 + gene_idx;
    for (; over_idx < genes_ptr->size()
        && (genes_ptr->get(over_idx).contig
          == genes_ptr->get(over_idx - 1).contig)
        && (genes_ptr->get(over_idx).coordinates.start
          <= genes_ptr->position_cummax()[over_idx - 1]); ++over_idx) { }
    // combined evidence for contig introns for these genes into a sorted vector
    // to process in sorted order
    auto stranded_evidence = OverGeneEvidence(
        exons, introns, gene_idx, over_idx, stranded);
    // loop over each strand/evidence for introns in each strand, add to result
    KnownContig contig = genes_ptr->get(gene_idx).contig;
    size_t prev_size = result_vec.size();
    for (const auto& [strand, evidence] : stranded_evidence) {
      AddIntronsFromEvidence(contig, strand, evidence, result_vec);
    }
    if (stranded_evidence.size() > 1) {
      // may not be in sorted order, so sort what we just added
      std::sort(result_vec.begin() + prev_size, result_vec.end());
    }
    // update gene_idx to over_idx
    gene_idx = over_idx;
  }
  // get final result
  return ContigIntrons{std::move(result_vec)};
}

}  // namespace majiq
