/**
 * TranscriptModels.cpp
 *
 * Implementation of ToSpliceGraph() for TranscriptModels class
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "TranscriptModels.hpp"

#include <utility>
#include <algorithm>
#include <sstream>
#include <stdexcept>


namespace majiq {
// create splicegraph with provided majiq gene models
SpliceGraph TranscriptModels::ToSpliceGraph(bool process_ir) const {
  // build exons, junctions, introns per gene
  std::vector<Exon> exons;
  std::vector<GeneJunction> junctions;
  std::vector<GeneIntron> introns;
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
          << gene.get().gene_id();
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
        << gene.get().gene_id();
      throw std::logic_error(oss.str());
    }
  }  // end iteration over all genes and their transcripts
  return SpliceGraph{contigs_, genes_,
    std::make_shared<Exons>(std::move(exons)),
    std::make_shared<GeneJunctions>(std::move(junctions)),
    std::make_shared<GeneIntrons>(std::move(introns))};
}
}  // namespace majiq
