/**
 * SpliceGraph.hpp
 *
 * Track different splicegraph elements
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_SPLICEGRAPH_HPP
#define MAJIQ_SPLICEGRAPH_HPP

#include <iostream>
#include <set>
#include <utility>
#include <vector>
#include <memory>
#include <algorithm>

#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "GeneJunctions.hpp"
#include "Exons.hpp"
#include "GeneIntrons.hpp"
#include "MajiqConstants.hpp"
#include "PassedJunctions.hpp"
#include "ExonConnections.hpp"


namespace majiq {

namespace detail {
/**
 * Infer exons from passed iterators for a single gene, add to dst
 *
 * @param exons_begin, exons_end range of exons to infer from
 * @param junctions_begin, junctions_end range of junctions to infer from
 * @param dst where inferred exons will be inserted (push_back)
 */
inline void GeneInferExons(
    Exons::const_iterator exons_begin,
    Exons::const_iterator exons_end,
    GeneJunctions::const_iterator junctions_begin,
    GeneJunctions::const_iterator junctions_end,
    std::vector<Exon>& dst);

}  // namespace detail

class SpliceGraph {
 protected:
  std::shared_ptr<Contigs> contigs_;
  std::shared_ptr<Genes> genes_;
  std::shared_ptr<Exons> exons_;
  std::shared_ptr<GeneJunctions> junctions_;
  std::shared_ptr<GeneIntrons> introns_;
  std::shared_ptr<ExonConnections> exon_connections_;

 public:
  // access non const pointers for use by pybind11 interface...
  std::shared_ptr<Contigs> contigs() { return contigs_; }
  std::shared_ptr<Genes> genes() { return genes_; }

  std::shared_ptr<Exons> exons() { return exons_; }

  std::shared_ptr<GeneJunctions> junctions() { return junctions_; }

  std::shared_ptr<GeneIntrons> introns() { return introns_; }

  std::shared_ptr<ExonConnections> exon_connections() {
    return exon_connections_;
  }

 private:
  template <typename ConnectionsT>
  static std::shared_ptr<ConnectionsT> ConnectedToExons(
      const std::shared_ptr<ConnectionsT>& connections,
      const std::shared_ptr<Exons>& exons) {
    connections->connect_exons(exons);
    return connections;
  }

 public:
  // constructors
  SpliceGraph(
      const std::shared_ptr<Contigs>& contigs,
      const std::shared_ptr<Genes>& genes,
      const std::shared_ptr<Exons>& exons,
      const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<GeneIntrons>& introns)
      : contigs_{contigs},
        genes_{genes},
        exons_{exons},
        junctions_{ConnectedToExons(junctions, exons)},
        introns_{ConnectedToExons(introns, exons)},
        exon_connections_{
          std::make_shared<ExonConnections>(exons_, introns_, junctions_)} { }
  SpliceGraph(const SpliceGraph& sg) = default;
  SpliceGraph(SpliceGraph&& sg) = default;
  SpliceGraph& operator=(const SpliceGraph& sg) = default;
  SpliceGraph& operator=(SpliceGraph&& sg) = default;

  static Exons InferExons(const Exons& source, const GeneJunctions& junctions) {
    if (source.parents() != junctions.parents()) {
      throw std::invalid_argument(
          "InferExons input exons and junctions do not share genes");
    }
    std::vector<Exon> result_vec;
    result_vec.reserve(source.size());  // will have at least same number exons
    for (size_t gene_idx = 0; gene_idx < source.parents()->size(); ++gene_idx) {
      // add exons inferred for each gene
      detail::GeneInferExons(
          source.begin_parent(gene_idx), source.end_parent(gene_idx),
          junctions.begin_parent(gene_idx), junctions.end_parent(gene_idx),
          result_vec);
    }
    return Exons{source.parents(), std::move(result_vec)};
  }

  GroupJunctionsGenerator MakeGroupGenerator() {
    return GroupJunctionsGenerator(junctions_, exons_);
  }
  PassedJunctionsGenerator MakePassedGenerator() {
    return PassedJunctionsGenerator(junctions_);
  }

  SpliceGraph Copy() {
    auto copy_junctions = std::make_shared<GeneJunctions>(*junctions_);
    auto copy_introns = std::make_shared<GeneIntrons>(*introns_);
    // NOTE: we don't have to copy contigs, genes, exons, which are effectively
    // immutable (can add to contigs, but will not affect what matters to copy)
    return SpliceGraph{
      contigs_, genes_, exons_, copy_junctions, copy_introns};
  }

  // to be declared later
  friend inline bool operator==(
      const SpliceGraph& x, const SpliceGraph& y) noexcept;
  friend inline std::ostream& operator<<(
      std::ostream&, const SpliceGraph&) noexcept;
};

// equality of splicegraphs
inline bool operator==(const SpliceGraph& x, const SpliceGraph& y) noexcept {
  return x.contigs_ == y.contigs_
    && x.genes_ == y.genes_
    && *x.exons_ == *y.exons_
    && *x.junctions_ == *y.junctions_
    && *x.introns_ == *y.introns_;
}

// how to represent splicegraph in output stream
inline std::ostream& operator<<(
    std::ostream& os, const SpliceGraph& sg) noexcept {
  os << "SpliceGraph<"
    << sg.contigs_->size() << " contigs, "
    << sg.genes_->size() << " genes, "
    << sg.exons_->size() << " exons, "
    << sg.junctions_->size() << " junctions, "
    << sg.introns_->size() << " introns"
    << ">";
  return os;
}

// implementation of detail::GeneInferExons (infer exons for single gene)
void detail::GeneInferExons(
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
  ExonIntervalT coordinates, annotated;

  // how do we update this state?
  auto past_naked_acceptors = [&dst, &gene, &naked_acceptors](position_t x) {
    // if x more than MAX_DENOVO_DIFFERENCE, output naked_acceptors, return true
    bool result = naked_acceptors.empty()
      || (naked_acceptors.back() + MAX_DENOVO_DIFFERENCE < x);
    if (result) {
      // add the naked acceptors and clear
      std::transform(naked_acceptors.begin(), naked_acceptors.end(),
          std::back_inserter(dst), [&gene](position_t y) {
          return Exon{gene, ExonIntervalT{y, -1}, Exon::MakeDenovo{}}; });
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
              gene, ExonIntervalT{-1, position}, Exon::MakeDenovo{});
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

}  // namespace majiq

#endif  // MAJIQ_SPLICEGRAPH_HPP
