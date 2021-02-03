/**
 * GFF3.hpp
 *
 * Process GFF3 file for MAJIQ genes/transcripts/exons
 *
 * Copyright 2020 <University of Pennsylvania
 */
#ifndef MAJIQ_GFF3_HPP
#define MAJIQ_GFF3_HPP

#include <memory>
#include <set>
#include <map>
#include <unordered_map>
#include <string>
#include <utility>
#include <variant>
#include <optional>

#include "Contigs.hpp"
#include "Genes.hpp"
#include "TranscriptModels.hpp"

namespace majiq {
namespace gff3 {
// number of columns in GFF3
constexpr size_t NUM_COLUMNS = 9;
// columns of interest
constexpr size_t COL_SEQID = 0;
constexpr size_t COL_TYPE = 2;
constexpr size_t COL_START = 3;
constexpr size_t COL_END = 4;
constexpr size_t COL_STRAND = 6;
constexpr size_t COL_ATTRIBUTES = 8;


using feature_id_t = std::string;
using gene_or_ancestor_id_t = std::variant<size_t, feature_id_t>;
using skipped_features_ct_t = std::map<std::string, unsigned int>;

struct GFF3TranscriptModels {
  const TranscriptModels models_;
  const skipped_features_ct_t skipped_transcript_type_ct_;
  const skipped_features_ct_t skipped_gene_type_ct_;

  GFF3TranscriptModels(
      TranscriptModels&& models,
      skipped_features_ct_t&& skipped_transcript_type_ct,
      skipped_features_ct_t&& skipped_gene_type_ct)
      : models_{models},
        skipped_transcript_type_ct_{skipped_transcript_type_ct},
        skipped_gene_type_ct_{skipped_gene_type_ct} {
  }
  GFF3TranscriptModels() = delete;
  GFF3TranscriptModels(const GFF3TranscriptModels&) = default;
  GFF3TranscriptModels(GFF3TranscriptModels&&) = default;
  GFF3TranscriptModels& operator=(const GFF3TranscriptModels&) = delete;
  GFF3TranscriptModels& operator=(GFF3TranscriptModels&&) = delete;
};

class GFF3ExonHierarchy {
 private:
  // contigs and genes recognized so far
  const std::shared_ptr<Contigs> contigs_;
  const std::shared_ptr<Genes> genes_;
  // feature_genes_: feature id --> size_t if known gene, otherwise some
  // ancestor (or self if top-level feature). Updated toward top-level feature
  // or first gene ancestor using disjoint set find
  std::unordered_map<feature_id_t, gene_or_ancestor_id_t> feature_genes_;
  // the set of exons that are children to a given feature id
  std::unordered_map<feature_id_t, transcript_exons_t> feature_exons_;
  // track feature types that aren't exons. Empty if accepted transcript/gene
  std::unordered_map<feature_id_t, std::optional<std::string>>
    feature_types_;

  /**
   * clear out non-const values when done using it (i.e. after std::move)
   */
  void clear() {
    feature_genes_.clear();
    feature_exons_.clear();
    feature_types_.clear();
  }

  /**
   * Given id, return gene or ancestor, with indicator if ancestor is top level
   *
   * Return <top_level, gene_or_oldest_ancestor>: if gene, top level. If oldest
   * ancestor, top_level indicates if top level or just missing information
   * about oldest ancestor (unable to follow to end at this time)
   */
  std::pair<bool, gene_or_ancestor_id_t> GetFeatureGene(feature_id_t id);

 public:
  // convert gff3 exon hierarchy to majiq transcript exons for processing
  friend GFF3TranscriptModels ToTranscriptModels(GFF3ExonHierarchy&&);

  /**
   * Load GFF3ExonHierarchy from specified input path
   */
  explicit GFF3ExonHierarchy(const std::string& gff3_filename);
  // default or deleted constructors/operators
  GFF3ExonHierarchy() = delete;
  GFF3ExonHierarchy(const GFF3ExonHierarchy&) = default;
  GFF3ExonHierarchy(GFF3ExonHierarchy&&) = default;
  GFF3ExonHierarchy& operator=(const GFF3ExonHierarchy&) = delete;
  GFF3ExonHierarchy& operator=(GFF3ExonHierarchy&&) = delete;
};

GFF3TranscriptModels ToTranscriptModels(GFF3ExonHierarchy&&);


}  // namespace gff3
}  // namespace majiq


#endif  // MAJIQ_GFF3_HPP
