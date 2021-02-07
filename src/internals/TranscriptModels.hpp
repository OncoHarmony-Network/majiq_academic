/**
 * TranscriptModels.hpp
 *
 * 3-level hierarchy of genes, which have transcripts, which have exons, and
 * how to convert them into a splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_TRANSCRIPTMODELS_HPP
#define MAJIQ_TRANSCRIPTMODELS_HPP

#include <memory>
#include <vector>
#include <set>

#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "SpliceGraph.hpp"


namespace majiq {
using transcript_exons_t = std::set<ClosedInterval>;  // in sorted order

class TranscriptModels {
 public:
  const std::shared_ptr<Contigs> contigs_;
  const std::shared_ptr<Genes> genes_;
  const std::vector<std::vector<transcript_exons_t>> gene_transcript_exons_;

  size_t size() const noexcept { return gene_transcript_exons_.size(); }
  SpliceGraph ToSpliceGraph(bool process_ir) const;

  // constructors
  TranscriptModels() = delete;
  TranscriptModels(const TranscriptModels&) = default;
  TranscriptModels(TranscriptModels&&) = default;
  TranscriptModels& operator=(const TranscriptModels&) = delete;
  TranscriptModels& operator=(TranscriptModels&&) = delete;
  TranscriptModels(
      const std::shared_ptr<Contigs>& contigs,
      const std::shared_ptr<Genes>& genes,
      std::vector<std::vector<transcript_exons_t>>&& gene_transcript_exons)
      : contigs_{contigs},
        genes_{genes},
        gene_transcript_exons_{gene_transcript_exons} {
    if (contigs_ == nullptr) {
      throw std::runtime_error("TranscriptModels needs non-null contigs");
    } else if (genes_ == nullptr) {
      throw std::runtime_error("TranscriptModels needs non-null genes");
    } else if (genes_->size() != gene_transcript_exons_.size()) {
      throw std::runtime_error("size mismatch in TranscriptModels");
    }
  }
};
}  // namespace majiq

#endif  // MAJIQ_TRANSCRIPTMODELS_HPP
