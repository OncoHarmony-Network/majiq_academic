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
#include <vector>
#include <memory>
#include <algorithm>

#include "Interval.hpp"
#include "Contigs.hpp"
#include "Genes.hpp"
#include "Junctions.hpp"
#include "Exons.hpp"
#include "Introns.hpp"


namespace majiq {
class SpliceGraph {
 protected:
  std::shared_ptr<Contigs> contigs_;
  std::shared_ptr<Genes> genes_;
  std::shared_ptr<Exons> exons_;
  std::shared_ptr<GeneJunctions> junctions_;
  std::shared_ptr<Introns> introns_;
  // NOTE: if we add additional objects with KnownGene, update sort()

 public:
  // add to splicegraph
  void AddGene(const geneid_t& geneid, const seqid_t& seqid,
      const position_t& start, const position_t& end, bool strand_forward,
      genename_t genename) {
    genes_->add(Gene{
        contigs_->make_known(seqid), ClosedInterval{start, end},
        strand_forward ? GeneStrandness::FORWARD : GeneStrandness::REVERSE,
        geneid, genename});
    return;
  }

  // access non const pointers for use by pybind11 interface...
  std::shared_ptr<Contigs> contigs() { return contigs_; }
  std::shared_ptr<Genes> genes() { return genes_; }
  std::shared_ptr<Exons> exons() { return exons_; }
  std::shared_ptr<GeneJunctions> junctions() { return junctions_; }
  std::shared_ptr<Introns> introns() { return introns_; }

  // constructors
  SpliceGraph(
      std::shared_ptr<Contigs>&& contigs,
      std::shared_ptr<Genes>&& genes,
      std::shared_ptr<Exons>&& exons,
      std::shared_ptr<GeneJunctions>&& junctions,
      std::shared_ptr<Introns>&& introns)
      : contigs_{contigs},
        genes_{genes},
        exons_{exons},
        junctions_{junctions},
        introns_{introns} {
  }
  SpliceGraph(const SpliceGraph& sg) = default;
  SpliceGraph(SpliceGraph&& sg) = default;
  SpliceGraph& operator=(const SpliceGraph& sg) = default;
  SpliceGraph& operator=(SpliceGraph&& sg) = default;

  // to be declared later
  friend bool operator==(const SpliceGraph& x, const SpliceGraph& y) noexcept;
  friend std::ostream& operator<<(std::ostream&, const SpliceGraph&) noexcept;
};

// equality of splicegraphs
bool operator==(const SpliceGraph& x, const SpliceGraph& y) noexcept {
  return x.contigs_ == y.contigs_
    && x.genes_ == y.genes_
    && *x.exons_ == *y.exons_
    && *x.junctions_ == *y.junctions_
    && *x.introns_ == *y.introns_;
}

// how to represent splicegraph in output stream
std::ostream& operator<<(std::ostream& os, const SpliceGraph& sg) noexcept {
  os << "SpliceGraph<"
    << sg.contigs_->size() << " contigs, "
    << sg.genes_->size() << " genes, "
    << sg.exons_->size() << " exons, "
    << sg.junctions_->size() << " junctions, "
    << sg.introns_->size() << " introns"
    << ">";
  return os;
}
}  // namespace majiq

#endif  // MAJIQ_SPLICEGRAPH_HPP
