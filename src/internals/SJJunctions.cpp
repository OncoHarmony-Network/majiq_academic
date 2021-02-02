/**
 * SJJunctions.cpp
 *
 * Implementation of nontrivial functions in SJJunctions.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "SJJunctions.hpp"

namespace majiq {

void GroupJunctionsGenerator::AddExperiment(const SJJunctions& experiment) {
  for (const SJJunction& junction : experiment) {
    // need min positions to pass no matter what
    if (junction.data.numpos < minpos_) { continue; }
    // otherwise passed in general vs for denovo
    const bool minreads_passed = junction.data.numreads >= minreads_;
    const bool mindenovo_passed = junction.data.numreads >= mindenovo_;
    if (mindenovo_passed || minreads_passed) {
      // need to update counts_
      size_t contig_idx = contigs_->add(junction.contig.get());
      if (contig_idx >= counts_.size()) {
        // make sure index available for contig
        counts_.resize(contig_idx + 1);
      }
      // get reference to value for contig/coordinates/strand to add
      // make two copies of unstranded junctions to handle experiments with
      // mixed strandedness (unstranded/stranded)
      std::vector<GeneStrandness> strands;
      if (junction.strand == GeneStrandness::AMBIGUOUS) {
        // both strands in sorted order
        if constexpr(GeneStrandness::FORWARD < GeneStrandness::REVERSE) {
          strands = {GeneStrandness::FORWARD, GeneStrandness::REVERSE};
        } else {
          strands = {GeneStrandness::REVERSE, GeneStrandness::FORWARD};
        }
      } else {
        strands = {junction.strand};
      }
      for (GeneStrandness strand : strands) {
        auto& counts = counts_[contig_idx][
          std::make_pair(junction.coordinates, strand)];
        if (minreads_passed) {
          ++(counts.numpassed);
        }
        if (mindenovo_passed) {
          ++(counts.numdenovo);
        }
      }
    }  // done adding junction that passed
  }  // end loop over junctions
  // update num_experiments_
  ++num_experiments_;
}
}  // namespace majiq
