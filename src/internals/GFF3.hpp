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

using feature_id_t = std::string;
using skipped_features_ct_t = std::map<std::string, unsigned int>;
constexpr unsigned char _FLAG_IS_EXON = 1 << 0;
constexpr unsigned char _FLAG_IN_HIERARCHY = 1 << 1;
constexpr unsigned char _FLAG_ACCEPT_TRANSCRIPT = 1 << 2;
constexpr unsigned char _FLAG_SILENT = 1 << 3;
constexpr unsigned char _FLAG_IS_GENE = 1 << 4;
enum class FeatureType : unsigned char {
  EXON = _FLAG_IS_EXON,
  ACCEPT_GENE
    = _FLAG_IN_HIERARCHY | _FLAG_ACCEPT_TRANSCRIPT
    | _FLAG_IS_GENE | _FLAG_SILENT,
  ACCEPT_TRANSCRIPT = _FLAG_IN_HIERARCHY | _FLAG_ACCEPT_TRANSCRIPT,
  REJECT_SILENT = _FLAG_IN_HIERARCHY | _FLAG_SILENT,
  REJECT_OTHER = _FLAG_IN_HIERARCHY,
  HARD_SKIP = 0
};
inline bool type_in_hierarchy(FeatureType x) {
  return static_cast<unsigned char>(x) & _FLAG_IN_HIERARCHY;
}
inline bool type_accepted_transcript(FeatureType x) {
  return static_cast<unsigned char>(x) & _FLAG_ACCEPT_TRANSCRIPT;
}
inline bool type_is_silent(FeatureType x) {
  return static_cast<unsigned char>(x) & _FLAG_SILENT;
}
using featuretype_map_t = std::map<std::string, FeatureType>;
using featuretype_info_t = std::pair<FeatureType, std::string>;

// throw exception if featuretype_map_t has no chance of succeeding
inline void assert_featuretype_map_plausible(
    const featuretype_map_t& gff3_types) {
  bool no_gene = true;  // it needs to have at least one ACCEPT_GENE
  bool no_exon = true;  // it needs to have at least one EXON
  for (const auto& x : gff3_types) {
    if (x.second == FeatureType::ACCEPT_GENE) {
      no_gene = false;
    } else if (x.second == FeatureType::EXON) {
      no_exon = false;
    }
  }
  if (no_gene || no_exon) {
    throw std::invalid_argument(
        "GFF3 feature map needs at least one ACCEPT_GENE and one EXON");
  }
  return;
}


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
  std::shared_ptr<Genes> genes_;
  // feature_genes_: feature id --> size_t if known gene, otherwise some
  // ancestor (or self if top-level feature). Updated toward top-level feature
  // or first gene ancestor using disjoint set find
  std::unordered_map<feature_id_t, feature_id_t> feature_genes_;
  // the set of exons that are children to a given feature id
  std::unordered_map<feature_id_t, transcript_exons_t> feature_exons_;
  // track feature types that aren't exons. Track names when OTHER
  std::unordered_map<feature_id_t, featuretype_info_t> feature_types_;

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
  std::pair<bool, feature_id_t> GetFeatureGene(const feature_id_t& id);

 public:
  // convert gff3 exon hierarchy to majiq transcript exons for processing
  friend GFF3TranscriptModels ToTranscriptModels(GFF3ExonHierarchy&&);

  /**
   * Load GFF3ExonHierarchy from specified input path
   */
  GFF3ExonHierarchy(const std::string& gff3_filename, const featuretype_map_t&);
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
