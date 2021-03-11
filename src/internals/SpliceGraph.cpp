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

#include "MajiqConstants.hpp"


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
    sites.emplace(eit->annotated_coordinates().start, SpliceT::ANN_EXON_START);
    sites.emplace(eit->annotated_coordinates().end, SpliceT::ANN_EXON_END);
  }

  // iterate over sites to build our result
  // define state: potential half-acceptors/extension, current exon.
  std::vector<position_t> naked_acceptors;
  ClosedInterval coordinates, annotated;

  // how do we update this state?
  auto past_naked_acceptors = [&dst, &gene, &naked_acceptors](position_t x) {
    // if x more than MAX_DENOVO_DIFFERENCE, output naked_acceptors, return true
    bool result = naked_acceptors.empty()
      || (naked_acceptors.back() + MAX_DENOVO_DIFFERENCE < x);
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
      || (coordinates.end + MAX_DENOVO_DIFFERENCE < x);
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
  return Exons{source.parents(), std::move(result_vec)};
}

}  // namespace majiq
