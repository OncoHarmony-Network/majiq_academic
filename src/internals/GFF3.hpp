/**
 * GFF3.hpp
 *
 * Simple interface to MAJIQ GFF3 file
 *
 * Copyright 2020 <University of Pennsylvania
 */
#ifndef MAJIQ_GFF3_HPP
#define MAJIQ_GFF3_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <memory>
#include <regex>
#include <stdexcept>
#include <cstdlib>
#include <boost/algorithm/string.hpp>

#include "Contigs.hpp"
#include "Genes.hpp"
#include "GeneJunctions.hpp"
#include "Introns.hpp"
#include "Exons.hpp"
#include "SpliceGraph.hpp"

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

enum class Feature : char {
  GENE,
  TRANSCRIPT,
  EXON,
  OTHER
};


class SpliceGraphBuilder {
 private:
  // contigs and genes recognized
  std::shared_ptr<Contigs> contigs_;
  std::shared_ptr<Genes> genes_;
  // have to handle transcripts, which we flatten
  using transcript_id_t = std::string;
  using transcript_exons_t = std::set<ClosedInterval>;
  // transcripts map to known or unknown genes
  std::map<transcript_id_t, size_t> known_transcripts_;  // map to gene_idx
  std::map<transcript_id_t, geneid_t> orphan_transcripts_;  // map to unseen id
  // exons map to known gene_idx/transcript or otherwise
  std::vector<std::map<transcript_id_t, transcript_exons_t>> known_exons_;
  std::map<transcript_id_t, transcript_exons_t> orphan_exons_;

  // static conversion from COL_TYPE to feature
  static const std::map<std::string, Feature> type_feature_map;
  static Feature TypeToFeature(std::string x) {
    auto match = type_feature_map.find(x);
    return match == type_feature_map.end() ? Feature::OTHER : match->second;
  }
  // extract IDs from attribute strings
  static const std::vector<std::regex> regex_geneid_vec;
  static const std::vector<std::regex> regex_genename_vec;
  static const std::vector<std::regex> regex_transcript_id_vec;
  static const std::vector<std::regex> regex_parent_vec;
  /**
   * @param patterns vector of capturing regexes -- first match is result, try
   * first pattern, then next if no match, etc.
   */
  inline static std::string attribute_value(
      const std::string& attributes, const std::vector<std::regex>& patterns) {
    std::smatch match;
    for (auto& pattern : patterns) {
      if (std::regex_search(attributes, match, pattern)) {
        return match[1];
      }
    }
    // if we get here, we don't have a match...
    std::ostringstream oss;
    oss << "Unable to obtain match for provided GFF3 attributes: "
      << attributes;
    throw std::runtime_error(oss.str());
  }
  inline static ClosedInterval get_interval(
      const std::vector<std::string>& record) {
    position_t start = std::atol(record[COL_START].c_str());
    position_t end = std::atol(record[COL_END].c_str());
    return ClosedInterval{start, end};
  }
  inline static GeneStrandness convert_strand(const std::string& col_strand) {
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

  // process a gene record
  void ProcessGeneRecord(const std::vector<std::string>& record) {
    // extract relevant information
    KnownContig contig = contigs_->make_known(record[COL_SEQID]);
    ClosedInterval interval = get_interval(record);
    GeneStrandness strand = convert_strand(record[COL_STRAND]);
    geneid_t geneid = attribute_value(
        record[COL_ATTRIBUTES], regex_geneid_vec);
    geneid_t genename = attribute_value(
        record[COL_ATTRIBUTES], regex_genename_vec);
    // make gene
    genes_->add(Gene{contig, interval, strand, geneid, genename});
    return;
  }

  // process a transcript record
  void ProcessTranscriptRecord(const std::vector<std::string>& record) {
    transcript_id_t id = attribute_value(
        record[COL_ATTRIBUTES], regex_transcript_id_vec);
    geneid_t geneid = attribute_value(
        record[COL_ATTRIBUTES], regex_parent_vec);
    std::optional<size_t> gene_idx = genes_->safe_gene_idx(geneid);
    if (gene_idx.has_value()) {
      known_transcripts_[id] = *gene_idx;
    } else {
      orphan_transcripts_[id] = geneid;
    }
  }

  // process an exon record
  void ProcessExonRecord(const std::vector<std::string>& record) {
    transcript_id_t transcript_id = attribute_value(
        record[COL_ATTRIBUTES], regex_parent_vec);
    ClosedInterval interval = get_interval(record);
    // is this a known gene-transcript?
    auto find_gene_idx = known_transcripts_.find(transcript_id);
    if (find_gene_idx == known_transcripts_.end()) {
      // this is an orphan exon
      orphan_exons_[transcript_id].insert(interval);
    } else {
      // get gene_idx and make sure known exons has record for gene
      size_t gene_idx = find_gene_idx->second;
      if (gene_idx >= known_exons_.size()) {
        known_exons_.resize(gene_idx + 1);
      }
      // add it to appropriate set of exons
      known_exons_[gene_idx][transcript_id].insert(interval);
    }
  }

  // dispatch record to appropriate feature handler
  void ProcessRecord(const std::vector<std::string>& record) {
    switch (TypeToFeature(record[COL_TYPE])) {
      case Feature::GENE:
        ProcessGeneRecord(record);
        break;
      case Feature::TRANSCRIPT:
        ProcessTranscriptRecord(record);
        break;
      case Feature::EXON:
        ProcessExonRecord(record);
        break;
      default:
        // do nothing
        break;
    }
    return;
  }

  // process GFF3 file
  void ProcessGFF3(const std::string& filename) {
    std::ifstream in{filename};
    if (!in) {
      std::ostringstream oss;
      oss << "Failed to open " << filename << " in ProcessGFF3";
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
      // otherwise, we should process this record
      ProcessRecord(record);
    }
  }

  /**
   * Orphans might be known after processing full annotation.
   * Either add them to known transcripts or report them as missing.
   *
   * Clear known_transcripts_ once done as all transcripts have been mapped to
   * genes (or reported as orphans)
   */
  void ProcessOrphans() {
    // get iterators over orphans
    auto it_transcripts = orphan_transcripts_.begin();
    auto it_exons = orphan_exons_.begin();
    // what do we do when we are processing a transcript?
    auto process_transcript = [this, &it_transcripts]() {
      // do we recognize this gene now?
      std::optional<size_t> gene_idx
        = genes_->safe_gene_idx(it_transcripts->second);
      if (gene_idx.has_value()) {
        known_transcripts_[it_transcripts->first] = *gene_idx;
      } else {
        std::cerr << "Orphan transcript " << it_transcripts->first
          << " with unfound gene " << it_transcripts->second << "\n";
      }
      // remove from orphan_transcripts_
      orphan_transcripts_.erase(it_transcripts++);
    };
    auto process_exons = [this, &it_exons]() {
      // do we recognize the transcript now?
      auto known = known_transcripts_.find(it_exons->first);
      if (known == known_transcripts_.end()) {
        // the transcript is not known
        std::cerr << it_exons->second.size()
          << " orphan exons with unfound transcript "
          << it_exons->first << "\n";
      } else {
        // merge orphan exons with set of exons that were known
        known_exons_[known->second][it_exons->first].merge(it_exons->second);
      }
      // remove from orphan_exons_
      orphan_exons_.erase(it_exons++);
    };
    // can process in a merge-like fashion when we have both remaining
    while (it_transcripts != orphan_transcripts_.end()
        && it_exons != orphan_exons_.end()) {
      // transcripts have priority over exons if they match transcript_id
      if (it_transcripts->first <= it_exons->first) {
        process_transcript();  // updates it_transcripts
      } else {
        process_exons();  // updates it_exons
      }
    }
    // process any remaining transcripts/exons
    for (; it_transcripts != orphan_transcripts_.end();) {
      process_transcript();
    }
    for (; it_exons != orphan_exons_.end();) {
      process_exons();
    }
    // clear known transcripts
    known_transcripts_.clear();
    return;
  }

  // NOTE: right now we are ignoring transcript start/ends
  void FlattenGene(
      KnownGene gene,
      // flatten the gene models
      const std::map<transcript_id_t, transcript_exons_t>& gene_models,
      // into the following vectors
      std::vector<Exon>* exons,
      std::vector<GeneJunction>* junctions,
      std::vector<Intron>* introns,
      // whether we do anything with introns...
      bool process_ir) {
    if (gene_models.empty()) {
      return;
    }
    // build up a set of junctions for this gene so we can add in order
    std::set<OpenInterval> gene_junctions;
    // exons and introns sorted and unique by construction
    std::vector<ClosedInterval> gene_exons;
    std::vector<ClosedInterval> gene_introns;
    // array of splice sites (we count duplicates)
    // remember that start comes before end. So pair<position, is_end>.
    std::vector<std::pair<position_t, bool>> splice_sites;

    // collect junctions and splicesites (TODO(jaicher): tx start/ends too)
    for (const auto& [tx_id, tx_exons] : gene_models) {
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
    for (const auto [pos, is_end] : splice_sites) {
      // update how many starts we've had more than ends, check it's valid
      if ((nopen += is_end ? -1 : 1) < 0) {
        // we have had more ends than starts. This should never happen
        std::ostringstream oss;
        oss << "More exon ends before starts when flattening gene "
          << gene.get().geneid;
        throw std::logic_error(oss.str());
      }
      if (is_end) {
        cur_end = pos;  // sorting --> pos >= cur_end
      } else {
        if (cur_start == EMPTY) {
          // this the first start in the gene
          cur_start = pos;
        } else if (cur_end != EMPTY) {
          // this is a start that follows an end, creating an exon
          gene_exons.push_back(ClosedInterval{cur_start, cur_end});
          if (process_ir && nopen > 1) {
            // there is an exon besides the one that just opened that's still
            // running from newly added exon. Treated as an annotated intron.
            gene_introns.push_back(ClosedInterval{cur_end + 1, pos - 1});
          }
          // clear cur_end, set cur_start
          cur_end = EMPTY;
          cur_start = pos;
        }  // end cases where first start of an exon
      }  // end handling of splice sites at start vs end of exon
    }  // end iteration over all splice sites
    // add final exon
    gene_exons.push_back(ClosedInterval{cur_start, cur_end});
    if (nopen != 0) {
      // we have had more starts than ends. This should never happen.
      std::ostringstream oss;
      oss << "Found more exon starts than ends when flattening gene "
        << gene.get().geneid;
      throw std::logic_error(oss.str());
    }

    // update junctions/exons/introns now past possible exceptions
    std::for_each(gene_junctions.begin(), gene_junctions.end(),
        [gene, &junctions](OpenInterval iv) {
          junctions->push_back(GeneJunction{gene, iv});
        });
    std::for_each(gene_exons.begin(), gene_exons.end(),
        [gene, &exons](ClosedInterval iv) {
          exons->push_back(Exon{gene, iv, Exon::DefaultAnnotated{}});
        });
    std::for_each(gene_introns.begin(), gene_introns.end(),
        [gene, &introns](ClosedInterval iv) {
          introns->push_back(Intron{gene, iv});
        });
    return;
  }

  SpliceGraph BuildSpliceGraph(bool process_ir) {
    // shared pointers
    std::vector<Exon> exons;
    std::vector<GeneJunction> junctions;
    std::vector<Intron> introns;
    // get sorted genes
    std::shared_ptr<Genes> sorted_genes = genes_->sorted();
    // for each gene in order to make result sorted in order
    for (size_t gene_idx = 0; gene_idx < known_exons_.size(); ++gene_idx) {
      // flatten the sorted genes in order
      const KnownGene g = (*sorted_genes)[gene_idx];
      // get the corresponding gene model
      auto& model = known_exons_[
        sorted_genes == genes_
        ? gene_idx
        : genes_->get_gene_idx(sorted_genes->get(gene_idx))
      ];
      FlattenGene(g, model, &exons, &junctions, &introns, process_ir);
      // clear out model
      model.clear();
    }
    // construct splicegraph
    SpliceGraph result{std::move(contigs_), std::move(sorted_genes),
      std::make_shared<Exons>(exons),
      std::make_shared<GeneJunctions>(junctions),
      std::make_shared<Introns>(introns)
    };
    // reset builder
    contigs_ = std::make_shared<Contigs>();
    genes_ = std::make_shared<Genes>();
    known_exons_.clear();
    known_exons_.shrink_to_fit();
    // return result
    return result;
  }

 public:
  SpliceGraphBuilder()
      : contigs_{std::make_shared<Contigs>()},
        genes_{std::make_shared<Genes>()} {
  }
  SpliceGraphBuilder(const SpliceGraphBuilder& x) = default;
  SpliceGraphBuilder(SpliceGraphBuilder&& x) = default;
  SpliceGraphBuilder& operator=(const SpliceGraphBuilder& x) = default;
  SpliceGraphBuilder& operator=(SpliceGraphBuilder&& x) = default;

  SpliceGraph from_gff3(const std::string& gff3_filename, bool process_ir) {
    ProcessGFF3(gff3_filename);  // builds gene/transcript/exon model
    ProcessOrphans();  // clears orphan transcripts/exons, transcript-gene map
    return BuildSpliceGraph(process_ir);  // resets contigs/genes/known exons
  }
};

const std::map<std::string, Feature> SpliceGraphBuilder::type_feature_map = {
  // accepted transcripts
  {"mRNA", Feature::TRANSCRIPT},
  {"transcript", Feature::TRANSCRIPT},
  {"lnc_RNA", Feature::TRANSCRIPT},
  {"miRNA", Feature::TRANSCRIPT},
  {"ncRNA", Feature::TRANSCRIPT},
  {"rRNA", Feature::TRANSCRIPT},
  {"scRNA", Feature::TRANSCRIPT},
  {"snRNA", Feature::TRANSCRIPT},
  {"snoRNA", Feature::TRANSCRIPT},
  {"tRNA", Feature::TRANSCRIPT},
  {"pseudogenic_transcript", Feature::TRANSCRIPT},
  {"C_gene_segment", Feature::TRANSCRIPT},
  {"D_gene_segment", Feature::TRANSCRIPT},
  {"J_gene_segment", Feature::TRANSCRIPT},
  {"V_gene_segment", Feature::TRANSCRIPT},
  {"unconfirmed_transcript", Feature::TRANSCRIPT},
  {"three_prime_overlapping_ncrna", Feature::TRANSCRIPT},
  // accepted genes
  {"gene", Feature::GENE},
  {"ncRNA_gene", Feature::GENE},
  {"pseudogene", Feature::GENE},
  {"ncRNA_gene", Feature::GENE},
  {"bidirectional_promoter_lncRNA", Feature::GENE},
  // accepted exons
  {"exon", Feature::EXON}
};
const std::vector<std::regex> SpliceGraphBuilder::regex_geneid_vec = {
  std::regex{"(?:^|;)ID=([^;]+)"},
  std::regex{"(?:^|;)gene_id=([^;]+)"}
};
const std::vector<std::regex> SpliceGraphBuilder::regex_genename_vec = {
  std::regex{"(?:^|;)Name=([^;]+)"},
  std::regex{"(?:^|;)gene_name=([^;]+)"}
};
const std::vector<std::regex> SpliceGraphBuilder::regex_transcript_id_vec = {
  std::regex{"(?:^|;)ID=([^;]+)"}
};
const std::vector<std::regex> SpliceGraphBuilder::regex_parent_vec
= {
  std::regex{"(?:^|;)Parent=([^;]+)"}
};

}  // namespace gff3
}  // namespace majiq

#endif  // MAJIQ_GFF3_HPP
