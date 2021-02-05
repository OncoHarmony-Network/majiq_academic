/**
 * GFF3.cpp
 *
 * Implementation of GFF3.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "GFF3.hpp"

#include <utility>
#include <stdexcept>
#include <sstream>
#include <set>
#include <unordered_set>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <regex>

#include <zstr/zstr.hpp>
#include <boost/algorithm/string.hpp>

namespace majiq {
namespace gff3 {

/**
 * Gets gene/top-level feature corresponding to input feature if possible
 *
 * Updates features in path to gene/top-level feature using path compression
 *
 * Return <top_level, gene_or_oldest_ancestor>: if gene, top level. If oldest
 * ancestor, top_level indicates if top level or just missing information
 * about oldest ancestor (unable to follow to end at this time)
 */
std::pair<bool, gene_or_ancestor_id_t> GFF3ExonHierarchy::GetFeatureGene(
    feature_id_t id) {
auto cur_parent = feature_genes_.find(id);
if (cur_parent != feature_genes_.end()) {
  // get gene or parent id that was found
  const auto& gene_or_parent_id = cur_parent->second;
  // if it is a gene_idx, we are done
  if (auto gene_idx_ptr = std::get_if<size_t>(&gene_or_parent_id)) {
    return std::make_pair(true, gene_or_ancestor_id_t{*gene_idx_ptr});
  }
  // otherwise, we have a parent id
  const auto& parent_id = std::get<feature_id_t>(gene_or_parent_id);
  // if parent is self, this is top level!
  if (id == parent_id) {
    return std::make_pair(true, gene_or_ancestor_id_t{parent_id});
  } else {
    // get result
    auto result = GetFeatureGene(parent_id);
    // bind to values
    const auto& [is_top_level, final_gene_or_ancestor] = result;
    // update where id points to
    feature_genes_[id] = final_gene_or_ancestor;
    // pass result on down
    return result;
  }
}
// otherwise, no parent so not top level
return std::make_pair(false, gene_or_ancestor_id_t{id});
}

/**
 * Consume parsed GFF3 exon hierarchy, make majiq exon hierarchy
 */
GFF3TranscriptModels ToTranscriptModels(GFF3ExonHierarchy&& x) {
  // get genes in sorted order
  std::shared_ptr<Genes> genes = x.genes_->sorted();
  bool x_genes_sorted = genes == x.genes_;
  // initialize vector of transcript/gene models per gene
  std::vector<std::vector<transcript_exons_t>> gene_transcript_exons(
      genes->size());
  // initialize counts of skipped types
  std::map<std::string, unsigned int> skipped_transcript_type_ct;
  std::map<std::string, unsigned int> skipped_gene_type_ct;
  std::unordered_set<feature_id_t> skipped_genes;

  // loop over feature exons, assigning them appropriately
  // (could erase as we went, but probably premature optimization)
  for (const auto& [parent_id, exons] : x.feature_exons_) {
    {
      // did we see parent_id?
      const auto feature_types_it = x.feature_types_.find(parent_id);
      if (feature_types_it == x.feature_types_.end()) {
        // throw error if we don't have parent at all
        std::ostringstream oss;
        oss << "Found exons with parent id '" << parent_id << "' never defined";
        throw std::runtime_error(oss.str());
      }
      const auto& [parent_type, parent_type_str] = feature_types_it->second;
      // acceptable as a transcript?
      if (!type_accepted_transcript(parent_type)) {
        // no, so we skip these exons
        if (!type_is_silent(parent_type)) {
          // we also report the type of this one
          ++skipped_transcript_type_ct[parent_type_str];
        }
        continue;
      }
    }

    // is it or one of its ancestors a gene?
    const auto& [top_level, final_gene_or_ancestor]
      = x.GetFeatureGene(parent_id);
    if (!top_level) {
      // throw error if we don't have final gene or top-level ancestor...
      std::ostringstream oss;
      oss << "Found exons with undefined ancestor (ancestor ID: '"
        << std::get<feature_id_t>(final_gene_or_ancestor)
        << "')";
      throw std::runtime_error(oss.str());
    } else if (auto skipped_id_ptr
        = std::get_if<feature_id_t>(&final_gene_or_ancestor)) {
      // we are skipping this gene and not processing these
      skipped_genes.insert(*skipped_id_ptr);
    } else {
      // index into original genes built when parsing gff3
      const auto& x_gene_idx = std::get<size_t>(final_gene_or_ancestor);
      // index into sorted genes
      const size_t gene_idx
        = x_genes_sorted
        ? x_gene_idx
        : genes->get_gene_idx(x.genes_->get(x_gene_idx));
      // steal exons from x.feature_exons_ for resulting gene_transcript_exons
      gene_transcript_exons[gene_idx].push_back(std::move(exons));
    }
  }  // done copying over all exons that belong to genes/transcripts

  // count types of skipped genes
  for (const auto& gene_id : skipped_genes) {
    const auto& [gene_type, gene_type_str] = x.feature_types_[gene_id];
    if (!type_is_silent(gene_type)) {
      ++skipped_gene_type_ct[gene_type_str];
    }
  }
  // clear state of non-const elements of x (since passed as rvalue reference)
  x.clear();
  // construct and return result
  return GFF3TranscriptModels{
    TranscriptModels{x.contigs_, genes, std::move(gene_transcript_exons)},
    std::move(skipped_transcript_type_ct), std::move(skipped_gene_type_ct)};
}


static const std::vector<std::regex> regex_parent_vec = {
  std::regex{"(?:^|;)Parent=([^;]+)"}
};
static const std::vector<std::regex> regex_id_vec = {
  std::regex{"(?:^|;)ID=([^;]+)"}
};
static const std::vector<std::regex> regex_gene_name_vec = {
  std::regex{"(?:^|;)Name=([^;]+)"},
  std::regex{"(?:^|;)gene_name=([^;]+)"}
};

inline std::optional<std::string> attribute_value(
    const std::string& attributes, const std::vector<std::regex>& patterns) {
  std::smatch match;
  for (auto& pattern : patterns) {
    if (std::regex_search(attributes, match, pattern)) {
      return match[1];
    }
  }
  // if we get here, we don't have a match...
  return std::nullopt;
}

inline ClosedInterval get_coordinates(const std::vector<std::string>& record) {
  position_t start = std::atol(record[COL_START].c_str());
  position_t end = std::atol(record[COL_END].c_str());
  return ClosedInterval{start, end};
}

inline GeneStrandness convert_strand(const std::string& col_strand) {
  switch (col_strand[0]) {
    case '+':
      return GeneStrandness::FORWARD;
    case '-':
      return GeneStrandness::REVERSE;
    default:
      std::ostringstream oss;
      oss << "Detected non-standard gene strand " << col_strand;
      throw std::runtime_error(oss.str());
  }
}

// load GFF3 exon hierarchy from input file
GFF3ExonHierarchy::GFF3ExonHierarchy(
    const std::string& gff3_filename, const featuretype_map_t& gff3_types)
    : contigs_{std::make_shared<Contigs>()}, genes_{std::make_shared<Genes>()} {
  zstr::ifstream in{gff3_filename};
  if (!in) {
    std::ostringstream oss;
    oss << "Failed to open "<< gff3_filename
      << " for constructing GFF3ExonHierarchy";
    throw std::runtime_error(oss.str());
  }
  std::string cur_line;  // holds current GFF3 line being processed
  while (std::getline(in, cur_line)) {
    if (cur_line.empty() || cur_line[0] == '#') {
      // skip blank/commented lines
      continue;
    }
    // split on tabs
    using boost::algorithm::split;
    using boost::algorithm::is_any_of;
    std::vector<std::string> record;  // holds each column
    split(record, cur_line, is_any_of("\t"));
    if (record.size() != NUM_COLUMNS) {
      std::ostringstream oss;
      oss << "GFF3 has record with incorrect number of columns:\n"
        << cur_line;
      throw std::runtime_error(oss.str());
    }

    // parse exon vs parent features (especially genes)
    FeatureType record_type = get_feature_type(record[COL_TYPE], gff3_types);
    if (record_type == FeatureType::EXON) {
      // assign coordinates for the exon to its parent
      const auto record_parent_opt
        = attribute_value(record[COL_ATTRIBUTES], regex_parent_vec);
      if (!record_parent_opt.has_value()) {
        // we expect exon to always have parent transcript/gene
        std::ostringstream oss;
        oss << "GFF3 has exon record without defined parent:\n" << cur_line;
        throw std::runtime_error(oss.str());
      } else {
        feature_exons_[*record_parent_opt].insert(get_coordinates(record));
      }
    } else if (type_in_hierarchy(record_type)) {
      // does non-exon record have an id?
      const auto record_id_opt
        = attribute_value(record[COL_ATTRIBUTES], regex_id_vec);
      if (!record_id_opt.has_value()) { continue; }
      const std::string& record_id = *record_id_opt;

      // save type for this record (more detailed if not silent)
      feature_types_[record_id] = std::make_pair(record_type,
          type_is_silent(record_type) ? std::string{} : record[COL_TYPE]);

      // get parent or gene_idx for the record
      const auto record_parent_opt
        = attribute_value(record[COL_ATTRIBUTES], regex_parent_vec);
      gene_or_ancestor_id_t record_parent
        = record_parent_opt.value_or(record_id);
      // special processing for genes, replace record_parent with gene_idx
      if (record_type == FeatureType::ACCEPT_GENE) {
        // extract contig information
        KnownContig contig = contigs_->make_known(Contig{record[COL_SEQID]});
        // construct gene components
        ClosedInterval coordinates = get_coordinates(record);
        GeneStrandness strand = convert_strand(record[COL_STRAND]);
        geneid_t geneid = record_id;
        // if no gene name defined, use gene id as name
        genename_t genename
          = attribute_value(record[COL_ATTRIBUTES], regex_gene_name_vec)
          .value_or(geneid);
        // add gene to genes_, get gene_idx for it
        size_t gene_idx
          = genes_->add(Gene{contig, coordinates, strand, geneid, genename});
        // update record_parent with gene_idx
        record_parent = gene_idx;
      }
      // initialize disjoint sets data structure with parent
      feature_genes_[record_id] = record_parent;
    }  // end handling non-exon feature (or exon previous ifelse block)
  }  // end iteration over lines in GFF3 file
}

}  // namespace gff3
}  // namespace majiq
