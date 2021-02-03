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
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <regex>

#include <boost/algorithm/string.hpp>

namespace majiq {
namespace gff3 {
// TODO(jaicher): move this out of GFF3.cpp
// create splicegraph with provided majiq gene models
SpliceGraph MajiqTranscriptExons::ToSpliceGraph(bool process_ir) const {
  // build exons, junctions, introns per gene
  std::vector<Exon> exons;
  std::vector<GeneJunction> junctions;
  std::vector<Intron> introns;
  // for each gene (in sorted order)
  for (size_t gene_idx = 0; gene_idx < size(); ++gene_idx) {
    // get transcripts associated with this gene, next one if no transcripts
    const std::vector<transcript_exons_t>& gene_models
      = gene_transcript_exons_[gene_idx];
    if (gene_models.empty()) { continue; }
    // get gene for this index
    KnownGene gene = (*genes_)[gene_idx];

    // array of splice sites, counting duplicates
    // start comes before end, so pair<position, is_end>
    std::vector<std::pair<position_t, bool>> splice_sites;
    // collect junctions and splicesites (TODO(jaicher): tx starts/ends too)
    {  // scope for sorted set of unique junctions
      std::set<OpenInterval> gene_junctions;
      for (const auto& tx_exons : gene_models) {
        position_t last_end = -2;
        for (const auto& iv : tx_exons) {
          // start and end of exon
          splice_sites.push_back({iv.start, false});
          splice_sites.push_back({iv.end, true});
          // junction from last exon
          // TODO(jaicher): handle transcript start/ends too
          if (last_end >= 0) {
            gene_junctions.insert(OpenInterval{last_end, iv.start});
          }
          last_end = iv.end;
        }
      }  // end loop over transcript exons
      // fill junctions with unique junctions for this gene
      std::for_each(gene_junctions.begin(), gene_junctions.end(),
          [gene, &junctions](OpenInterval iv) {
            junctions.emplace_back(gene, iv);
          });
    }

    // use splice sites to determine exons/introns
    std::sort(splice_sites.begin(), splice_sites.end());
    // Exons are split whenever a an exon end is followed by an exon start
    // An intron occurs whenever we haven't closed all starts with ends before
    // this happens.
    int nopen = 0;
    constexpr position_t EMPTY = -1;
    position_t cur_start = EMPTY;  // current putative exon start
    position_t cur_end = EMPTY;  // current putative exon end
    // each time that we add exon, we reset cur_end to EMPTY
    for (const auto& [pos, is_end] : splice_sites) {
      // update how many starts we've had more than ends.
      nopen += is_end ? -1 : 1;

      // if we have more ends than starts, something went wrong.
      if (nopen < 0) {
        std::ostringstream oss;
        oss << "More exon ends before starts when flattening gene "
          << gene.get().geneid;
        throw std::logic_error(oss.str());
      }

      // handle splice site end vs start
      if (is_end) {
        // current exon must go to at least current position because of sorting
        cur_end = pos;
      } else {
        if (cur_start == EMPTY) {
          // this the first start in the gene
          cur_start = pos;
        } else if (cur_end != EMPTY) {
          // this is a start that follows an end, creating an exon
          exons.emplace_back(
              gene, ClosedInterval{cur_start, cur_end},
              Exon::DefaultAnnotated{});
          if (process_ir && nopen > 1) {
            // there is an exon besides the one that just opened that's still
            // running from newly added exon. Treated as an annotated intron.
            introns.emplace_back(gene, ClosedInterval{cur_end + 1, pos - 1});
          }
          // clear cur_end, set cur_start
          cur_end = EMPTY;
          cur_start = pos;
        }  // end cases where first start of an exon
      }  // end handling of splice sites at start vs end of exon
    }  // end iteration over all splice sites

    // add final exon
    exons.emplace_back(
        gene, ClosedInterval{cur_start, cur_end}, Exon::DefaultAnnotated{});

    // check that we've closed all splice sites -- should have nopen = 0 at end
    if (nopen != 0) {
      // we have had more starts than ends. This should never happen.
      std::ostringstream oss;
      oss << "Found more exon starts than ends when flattening gene "
        << gene.get().geneid;
      throw std::logic_error(oss.str());
    }
  }  // end iteration over all genes and their transcripts
  return SpliceGraph{contigs_, genes_,
    std::make_shared<Exons>(exons),
    std::make_shared<GeneJunctions>(junctions),
    std::make_shared<Introns>(introns)};
}

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
MajiqTranscriptExons ToMajiqTranscriptExons(GFF3ExonHierarchy&& x) {
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
    // did we see parent_id? was it accepted (empty) or something else?
    const auto feature_types_it = x.feature_types_.find(parent_id);
    if (feature_types_it == x.feature_types_.end()) {
      // throw error if we don't have parent at all
      std::ostringstream oss;
      oss << "Found exons with parent id '" << parent_id << "' never defined";
      throw std::runtime_error(oss.str());
    } else if (feature_types_it->second.has_value()) {
      // skip this transcript, count unrecognized id towards skipped transcripts
      ++skipped_transcript_type_ct[*(feature_types_it->second)];
      continue;
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
    ++skipped_gene_type_ct[*(x.feature_types_[gene_id])];
  }
  // clear state of non-const elements of x (since passed as rvalue reference)
  x.clear();
  // construct and return result
  return MajiqTranscriptExons{x.contigs_, genes, gene_transcript_exons,
    skipped_transcript_type_ct, skipped_gene_type_ct};
}

// is the type an acceptable exon
inline bool is_exon(const std::string& type_str) {
  return type_str == "exon";
}

// types for acceptable genes
static const std::set<std::string> accepted_gene_types = {
  "gene",
  "ncRNA_gene",
  "pseudogene",
  "ncRNA_gene",
  "bidirectional_promoter_lncRNA"
};
inline bool is_acceptable_gene(const std::string& type_str) {
  return accepted_gene_types.count(type_str) > 0;
}
// types for acceptable transcripts (note: genes are considered acceptable
// transcripts too due to observations in RefSeq)
static const std::set<std::string> accepted_transcript_types = {
  "mRNA",
  "transcript",
  "lnc_RNA",
  "miRNA",
  "ncRNA",
  "rRNA",
  "scRNA",
  "snRNA",
  "snoRNA",
  "tRNA",
  "pseudogenic_transcript",
  "C_gene_segment",
  "D_gene_segment",
  "J_gene_segment",
  "V_gene_segment",
  "unconfirmed_transcript",
  "three_prime_overlapping_ncrna"
};
inline bool is_acceptable_transcript(
    bool is_gene, const std::string& type_str) {
  // genes are acceptable transcripts, otherwise search set
  return is_gene || (accepted_transcript_types.count(type_str) > 0);
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

inline ClosedInterval get_interval(const std::vector<std::string>& record) {
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
GFF3ExonHierarchy::GFF3ExonHierarchy(const std::string& gff3_filename)
    : contigs_{std::make_shared<Contigs>()}, genes_{std::make_shared<Genes>()} {
  std::ifstream in{gff3_filename};
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
    if (is_exon(record[COL_TYPE])) {
      const auto record_parent_opt
        = attribute_value(record[COL_ATTRIBUTES], regex_parent_vec);
      if (!record_parent_opt.has_value()) {
        // we expect exon to always have parent transcript/gene
        std::ostringstream oss;
        oss << "GFF3 has exon record without defined parent:\n" << cur_line;
        throw std::runtime_error(oss.str());
      } else {
        feature_exons_[*record_parent_opt].insert(get_interval(record));
      }
    } else {
      // does non-exon record have an id?
      const auto record_id_opt
        = attribute_value(record[COL_ATTRIBUTES], regex_id_vec);
      if (record_id_opt.has_value()) {
        const std::string& record_id = *record_id_opt;

        // save type for this record (empty for genes/transcripts)
        const bool is_gene = is_acceptable_gene(record[COL_TYPE]);
        const bool is_transcript = is_acceptable_transcript(
            is_gene, record[COL_TYPE]);
        feature_types_[record_id]
          = is_transcript
          ? std::optional<std::string>{} : std::make_optional(record[COL_TYPE]);

        // get parent or gene_idx for the record
        const auto record_parent_opt
          = attribute_value(record[COL_ATTRIBUTES], regex_parent_vec);
        gene_or_ancestor_id_t record_parent
          = record_parent_opt.value_or(record_id);
        // special processing for genes, replace record_parent with gene_idx
        if (is_gene) {
          // extract contig information
          KnownContig contig = contigs_->make_known(record[COL_SEQID]);
          // construct gene components
          ClosedInterval interval = get_interval(record);
          GeneStrandness strand = convert_strand(record[COL_STRAND]);
          geneid_t geneid = record_id;
          // if no gene name defined, use gene id as name
          genename_t genename
            = attribute_value(record[COL_ATTRIBUTES], regex_gene_name_vec)
            .value_or(geneid);
          // add gene to genes_, get gene_idx for it
          size_t gene_idx
            = genes_->add(Gene{contig, interval, strand, geneid, genename});
          // update record_parent with gene_idx
          record_parent = gene_idx;
        }
        // initialize disjoint sets data structure with parent
        feature_genes_[record_id] = record_parent;
      }  // end handling non-exon feature with ID
    }  // end handling non-exon feature (or exon previous ifelse block)
  }  // end iteration over lines in GFF3 file
}

}  // namespace gff3
}  // namespace majiq
