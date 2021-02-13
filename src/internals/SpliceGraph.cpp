/**
 * SpliceGraph.cpp
 *
 * Implementation of SpliceGraph functions
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "SpliceGraph.hpp"

#include <vector>
#include <set>
#include <utility>
#include <algorithm>


namespace majiq {

/**
 * Infer exons from passed iterators for a single gene, add to dst
 *
 * @param exons_begin, exons_end range of exons to infer from
 * @param junctions_begin, junctions_end range of junctions to infer from
 * @param dst where inferred exons will be inserted (push_back)
 */
void GeneInferExons(
    Exons::const_iterator exons_begin,
    Exons::const_iterator exons_end,
    GeneJunctions::const_iterator junctions_begin,
    GeneJunctions::const_iterator junctions_end,
    std::vector<Exon>& dst) {
  // special cases
  if (exons_begin == exons_end) {
    // no exons input --> no exons output
    return;
  } else if (junctions_begin == junctions_end) {
    // no junctions means exons cannot change, so just copy them
    std::copy(exons_begin, exons_end, std::back_inserter(dst));
    return;
  }

  // get gene all output exons will have
  const KnownGene& gene = exons_begin->gene;

  // reduce to inferring exons on using splice sites, but making special note
  // of annotated exons
  enum class SpliceT : unsigned char {
    ANN_EXON_START = 0,
    JUNCTION_END = 1,
    // we will ignore any junctions between start/end of annotated exons
    JUNCTION_START = 2,
    ANN_EXON_END = 3
  };
  using ss_t = std::pair<position_t, SpliceT>;
  // get sorted list of splicesites to work with
  std::set<ss_t> sites;
  for (auto jit = junctions_begin; jit != junctions_end; ++jit) {
    sites.emplace(jit->coordinates.start, SpliceT::JUNCTION_START);
    sites.emplace(jit->coordinates.end, SpliceT::JUNCTION_END);
  }
  for (auto eit = exons_begin; eit != exons_end; ++eit) {
    // only use full annotated exons
    if (eit->is_denovo() || !eit->is_full_exon()) { continue; }
    sites.emplace(eit->coordinates.start, SpliceT::ANN_EXON_START);
    sites.emplace(eit->coordinates.end, SpliceT::ANN_EXON_END);
  }

  // iterate over sites to build our result
  // define state: potential half-acceptors/extension, current exon.
  std::vector<position_t> naked_acceptors;
  ClosedInterval coordinates, annotated;

  // how do we update this state?
  auto past_naked_acceptors = [&dst, &gene, &naked_acceptors](position_t x) {
    // if x more than MAX_DENOVO_DIFFERENCE, output naked_acceptors, return true
    bool result = naked_acceptors.empty()
      || (naked_acceptors.back() + SpliceGraph::MAX_DENOVO_DIFFERENCE < x);
    if (result) {
      // add the naked acceptors and clear
      std::transform(naked_acceptors.begin(), naked_acceptors.end(),
          std::back_inserter(dst), [&gene](position_t y) {
          return Exon{gene, ClosedInterval{y, -1}, Exon::MakeDenovo{}}; });
      naked_acceptors.clear();
    }
    return result;
  };
  auto summarize_naked_acceptors = [&naked_acceptors]() -> position_t {
    // if naked acceptors but not past, we may only want the first value
    position_t result = naked_acceptors[0];
    naked_acceptors.clear();
    return result;
  };
  auto past_prev_exon = [&coordinates](position_t x) -> bool {
    return coordinates.is_invalid()
      || (coordinates.end + SpliceGraph::MAX_DENOVO_DIFFERENCE < x);
  };
  auto add_exon = [&gene, &coordinates, &annotated, &dst]() {
    if (!coordinates.is_invalid()) {
      dst.emplace_back(gene, coordinates, annotated);
      coordinates = {};
      annotated = {};
    }
    return;
  };

  // iterate over defined splicesites, updating state as appropriate
  for (auto it = sites.begin(); it != sites.end(); ++it) {
    const auto& [position, type] = *it;
    switch (type) {
      case SpliceT::ANN_EXON_START:
        add_exon();  // if we have an exon, add it
        // past_naked_acceptors true --> add as half exons, new exon start with
        // annotated boundary
        // past_naked_acceptors false --> exon extension backwards
        coordinates.start
          = past_naked_acceptors(position)
            ? position : summarize_naked_acceptors();
        annotated.start = position;
        // go to next annotated exon end, update end of coordinates, annotated
        it = std::find_if(it, sites.end(),
            [](const auto& x) { return x.second == SpliceT::ANN_EXON_END; });
        // annotated end serves as starting end for exon (we'll see if there is
        // extension from JUNCTION_END)
        coordinates.end = annotated.end = it->first;
        break;
      case SpliceT::JUNCTION_END:
        add_exon();  // if we have an exon, add it
        past_naked_acceptors(position);  // add past acceptors if too far away
        naked_acceptors.push_back(position);  // but this is new naked acceptor
        break;
      case SpliceT::JUNCTION_START:
        if (!past_prev_exon(position)) {
          // already have coordinates started, and this should extend it
          coordinates.end = position;
        } else if (!past_naked_acceptors(position)) {
          // denovo donor close enough to denovo acceptor for denovo exon
          coordinates.start = summarize_naked_acceptors();
          coordinates.end = position;
        } else {
          // we were past previous exon or acceptors
          add_exon();  // we add current exon (if it exists)
          dst.emplace_back(  // we add current value as well
              gene, ClosedInterval{-1, position}, Exon::MakeDenovo{});
        }
        break;
      case SpliceT::ANN_EXON_END:
        throw std::logic_error(
            "Found annotated exon end without preceding start");
        break;
    }
  }  // done looping over splice sites
  if (!naked_acceptors.empty()) {
    throw std::logic_error("InferExons is trying to lengthen a gene");
  }
  add_exon();  // add remaining exon, if any
  return;
}

Exons SpliceGraph::InferExons(
    const Exons& source, const GeneJunctions& junctions) {
  if (source.parents() != junctions.parents()) {
    throw std::invalid_argument(
        "InferExons input exons and junctions do not share genes");
  }
  std::vector<Exon> result_vec;
  result_vec.reserve(source.size());  // will have at least same number of exons
  for (size_t gene_idx = 0; gene_idx < source.parents()->size(); ++gene_idx) {
    // add exons inferred for each gene
    GeneInferExons(source.begin_parent(gene_idx), source.end_parent(gene_idx),
        junctions.begin_parent(gene_idx), junctions.end_parent(gene_idx),
        result_vec);
  }
  return Exons{std::move(result_vec)};
}

inline bool CloseToPrecedingFullExon(
    const Exons& exons, const KnownGene& gene, position_t x) {
  // get first exon past x
  auto it = exons.overlap_upper_bound(gene, x);
  --it;  // first exon behind that (overlaps or behind)
  // keep going backwards until full exon, too far away, or different gene
  for (; it != exons.begin(); --it) {
    if (it->gene != gene
        || it->is_full_exon()
        || (it->coordinates.last_pos()
          < x - SpliceGraph::MAX_DENOVO_DIFFERENCE)) {
      break;
    }
  }
  return (it->gene == gene
      && it->is_full_exon()
      && x >= it->coordinates.start
      && x <= it->coordinates.end + SpliceGraph::MAX_DENOVO_DIFFERENCE);
}
inline bool CloseToFollowingFullExon(
    const Exons& exons, const KnownGene& gene, position_t x) {
  // get first exon overlapping or past x
  auto it = exons.overlap_lower_bound(gene, x);
  // keep going forwards until full exon, too far away, or different gene
  for (; it != exons.end(); ++it) {
    if (it->gene != gene
        || it->is_full_exon()
        || (it->coordinates.first_pos()
          > x + SpliceGraph::MAX_DENOVO_DIFFERENCE)) {
      break;
    }
  }
  return (it != exons.end()
      && it->gene == gene
      && it->is_full_exon()
      && x >= it->coordinates.start - SpliceGraph::MAX_DENOVO_DIFFERENCE
      && x <= it->coordinates.end);
}

/**
 * Assign denovo junction pj to potential matched genes in range [first, last)
 *
 * Track genes in og that pj could be assigned to. So, we prioritize junctions
 * that are close (or overlapping) an exon in the gene, preferring genes where
 * the junction is close on both ends over just one end vs neither
 */
void AssignDenovoJunction(const PassedJunction& pj, const KnownGene& first,
    const KnownGene& last, const Exons& exons, std::vector<GeneJunction>& dst) {
  // when junction is close to exons on both junction ends as appropriate
  std::vector<KnownGene> close_x2;
  // when only one end of junction close to a gene exon
  std::vector<KnownGene> close_x1;
  // when within gene boundaries but not close to exons at all
  std::vector<KnownGene> far;

  // check all of the overgenes
  for (KnownGene gene = first; gene != last; ++gene) {
    if (pj.strand != gene.strand()
        || !IntervalSubsets(pj.coordinates, gene.coordinates())) {
      // we assume that contig matches (when selecting overgene)
      // but junction incompatible with individual genes on coordinates, strand
      continue;
    }
    // count how many of start/end are close to exon
    int n_close = 0;  // how many of start/end close to an exon?
    if (CloseToPrecedingFullExon(exons, gene, pj.coordinates.start)) {
      ++n_close;
    }
    if (CloseToFollowingFullExon(exons, gene, pj.coordinates.end)) {
      ++n_close;
    }
    // priority of this gene/where we track it
    auto& gene_priority_vec
      = n_close == 0 ? far : (n_close == 1 ? close_x1 : close_x2);
    gene_priority_vec.push_back(gene);
  }  // done looking at all potential genes
  // so get the highest priority vector of genes that is nonempty
  const auto& matched_genes
    = close_x2.empty() ? (close_x1.empty() ? far : close_x1) : close_x2;
  // add passed denovo junctions for these genes, if any
  std::transform(matched_genes.begin(), matched_genes.end(),
      std::back_inserter(dst),
      [&pj](const KnownGene& gene) {
      return GeneJunction{gene, pj.coordinates, true, true, false}; });
  return;
}

GeneJunctions SpliceGraph::AssignPassedJunctions(
    const GeneJunctions& source, const Exons& exons,
    const PassedJunctions& passed, bool add_denovo) {
  if (source.parents() != exons.parents()) {
    throw std::invalid_argument(
        "AssignPassedJunctions gene junctions/exons do not share genes");
  } else if (source.parents()->contigs() != passed.parents()) {
    throw std::invalid_argument(
        "AssignPassedJunctions genes/passed junctions do not share contigs");
  } else if (passed.empty()) {
    return source;  // if no passed junctions, should not change source
  }

  // otherwise, accumulate junctions in temporary vector
  std::vector<GeneJunction> result_vec;

  // loop over junctions per overgene in contig order vs passed junctions
  size_t last_idx = result_vec.size();  // start of previous overgene in result
  auto pj_it = passed.begin();
  for (auto og = source.parents()->overgene_begin();
      og != source.parents()->overgene_end();
      ++og) {
    // if pj_it not on same contig, skip to the contig
    // (processed so that any skipped junctions cannot be contained by any
    // overgene)
    if (pj_it->contig < og.contig) {
      // if contigs do not match, skip to og contig (processed so skipped
      // junctions cannot be contained by any overgene
      pj_it = passed.begin_parent(og.contig);
    }
    // go to first passed junction for which start can be in overgene
    // (any ignored cannot be contained by any overgene)
    pj_it = std::find_if(pj_it, passed.end(),
        [&og](const auto& x) {
        return og.contig < x.contig
        || og.coordinates.start <= x.coordinates.start; });
    // function to tell if this iterator is past the overgene
    auto pj_past_og = [&og, &pj_it, &passed]() -> bool {
      return pj_it == passed.end()
        || og.contig != pj_it->contig
        || og.coordinates.end < pj_it->coordinates.start;
    };

    // get gene junctions for this overgene in contig junction order
    std::vector<GeneJunction> og_junctions(
        source.begin_parent(og.first()),
        source.begin_parent(og.last()));
    std::sort(og_junctions.begin(), og_junctions.end(),
        detail::CompareContigStranded<GeneJunction>{});
    auto gj_it = og_junctions.begin();  // iterator into these junctions

    // when we either have original junctions or passed junctions left
    // merge from og_junctions and passed in contig order
    while (gj_it != og_junctions.end() || !pj_past_og()) {
      if (gj_it != og_junctions.end() && !pj_past_og() && *gj_it == *pj_it) {
        for (; gj_it != og_junctions.end() && *gj_it == *pj_it; ++gj_it) {
          // gene junction was passed
          gj_it->passed_build() = true;
          result_vec.push_back(*gj_it);
        }
        ++pj_it;  // done over gj_it that match pj_it
      }
      for (; gj_it != og_junctions.end() && (pj_past_og() || *gj_it < *pj_it);
          ++gj_it) {
        // gene junction from og_junctions was not passed, add it as is
        result_vec.push_back(*gj_it);
      }
      for (; !pj_past_og() && (gj_it == og_junctions.end() || *pj_it < *gj_it);
          ++pj_it) {
        // passed junction is denovo -- can/should we add it?
        if (add_denovo
            && pj_it->passed_denovo()
            && IntervalSubsets(pj_it->coordinates, og.coordinates)) {
          AssignDenovoJunction(
              *pj_it, og.first(), og.last(), exons, result_vec);
        }
      }
    }  // loop while any junctions remain

    // sort newly added junctions for this overgene
    std::sort(result_vec.begin() + last_idx, result_vec.end());
    last_idx = result_vec.size();
  }  // loop over overgenes to add junctions

  return GeneJunctions{std::move(result_vec)};
}

}  // namespace majiq
