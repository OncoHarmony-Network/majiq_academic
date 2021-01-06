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
  std::shared_ptr<Exons<>> exons_;
  std::shared_ptr<GeneJunctions<>> junctions_;
  std::shared_ptr<Introns<>> introns_;
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
  void AddExon(geneid_t geneid, position_t start, position_t end) {
    exons_->insert(Exon{genes_->known(geneid), ClosedInterval{start, end}});
    return;
  }
  void AddJunction(geneid_t geneid, position_t start, position_t end) {
    junctions_->insert(GeneJunction{
        genes_->known(geneid), OpenInterval{start, end}});
    return;
  }
  void AddIntron(geneid_t geneid, position_t start, position_t end) {
    introns_->insert(Intron{
        genes_->known(geneid), ClosedInterval{start, end}});
    return;
  }

  /**
   * Sort splicegraph in place
   */
  void sort() {
    // get sorted genes
    std::shared_ptr<Genes> sorted_genes = genes_->sorted();
    if (sorted_genes != genes_) {
      genes_ = sorted_genes;
      exons_->remap_genes(sorted_genes);
      junctions_->remap_genes(sorted_genes);
      introns_->remap_genes(sorted_genes);
      // NOTE: update any other objects with KnownGene
    } else {
      exons_->sort();
      junctions_->sort();
      introns_->sort();
    }
    // done
    return;
  }

  // access non const pointers for use by pybind11 interface...
  // NOTE: these sort the containers
  std::shared_ptr<Contigs> contigs() { return contigs_; }
  std::shared_ptr<Genes> genes() { sort(); return genes_; }
  std::shared_ptr<Exons<>> exons() { sort(); return exons_; }
  std::shared_ptr<GeneJunctions<>> junctions() { sort(); return junctions_; }
  std::shared_ptr<Introns<>> introns() { sort(); return introns_; }

  // constructors
  SpliceGraph()
      : contigs_{std::make_shared<Contigs>()},
        genes_{std::make_shared<Genes>()},
        exons_{std::make_shared<Exons<>>()},
        junctions_{std::make_shared<GeneJunctions<>>()},
        introns_{std::make_shared<Introns<>>()} {
  }
  // copy shares contigs/genes but copies exons, junctions, introns
  SpliceGraph(const SpliceGraph& sg)
      : contigs_{sg.contigs_},
        genes_{sg.genes_},
        exons_{std::make_shared<Exons<>>(*sg.exons_)},
        junctions_{std::make_shared<GeneJunctions<>>(*sg.junctions_)},
        introns_{std::make_shared<Introns<>>(*sg.introns_)} {
  }

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
