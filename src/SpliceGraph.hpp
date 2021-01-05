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
  // std::vector<Exon> exons_;
  Exons exons_;
  std::vector<GeneJunction> junctions_;
  std::vector<Intron> introns_;
  // NOTE: if we add additional objects with KnownGene, update sort()

 public:
  KnownGene AddGene(const Contig& contig, const ClosedInterval& interval,
      GeneStrandness strand, const geneid_t& geneid,
      const genename_t& genename) {
    return genes_->make_known(
        Gene{contigs_->make_known(contig),
          interval,
          strand,
          geneid,
          genename});
  }
  void AddExon(geneid_t geneid, ClosedInterval interval) {
    exons_.insert(Exon{genes_->known(geneid), interval});
    return;
  }
  const GeneJunction& AddJunction(geneid_t geneid, ClosedInterval interval) {
    junctions_.push_back(GeneJunction{genes_->known(geneid), interval});
    return junctions_.back();
  }
  const Intron& addIntron(geneid_t geneid, ClosedInterval interval) {
    introns_.push_back(Intron{genes_->known(geneid), interval});
    return introns_.back();
  }

  /**
   * Sort splicegraph in place
   */
  void sort() {
    // get sorted genes
    std::shared_ptr<Genes> sorted_genes = genes_->sorted();
    if (sorted_genes != genes_) {
      // update genes_
      genes_ = sorted_genes;
      // update exons
      exons_.remap_genes(sorted_genes);
      // we need to update gene_idx for all objects with KnownGene
      for (auto& x : junctions_) {
        x.gene = x.gene.remapped(sorted_genes);
      }
      for (auto& x : introns_) {
        x.gene = x.gene.remapped(sorted_genes);
      }
      // NOTE: update any other objects with KnownGene
    }
    // sort the vectors (NOTE update with additional vectors to sort)
    std::sort(junctions_.begin(), junctions_.end());
    std::sort(introns_.begin(), introns_.end());
    // done
    return;
  }

  // constructors
  SpliceGraph()
      : contigs_{std::make_shared<Contigs>()},
        genes_{std::make_shared<Genes>()},
        exons_{},
        junctions_{},
        introns_{} {
  }
  // copy shares contigs/genes but copies exons, junctions, introns
  SpliceGraph(const SpliceGraph& sg)
      : contigs_{sg.contigs_},
        genes_{sg.genes_},
        exons_{sg.exons_},
        junctions_{sg.junctions_.begin(), sg.junctions_.end()},
        introns_{sg.introns_.begin(), sg.introns_.end()} {
  }

  // to be declared later
  friend bool operator==(const SpliceGraph& x, const SpliceGraph& y) noexcept;
  friend std::ostream& operator<<(std::ostream&, const SpliceGraph&) noexcept;
};

// equality of splicegraphs
bool operator==(const SpliceGraph& x, const SpliceGraph& y) noexcept {
  return x.contigs_ == y.contigs_
    && x.genes_ == y.genes_
    && x.exons_ == y.exons_
    && x.junctions_ == y.junctions_
    && x.introns_ == y.introns_;
}

// how to represent splicegraph in output stream
std::ostream& operator<<(std::ostream& os, const SpliceGraph& sg) noexcept {
  os << "SpliceGraph<"
    << sg.contigs_->size() << " contigs, "
    << sg.genes_->size() << " genes, "
    << sg.exons_.size() << " exons, "
    << sg.junctions_.size() << " junctions, "
    << sg.introns_.size() << " introns"
    << ">";
  return os;
}
}  // namespace majiq

#endif  // MAJIQ_SPLICEGRAPH_HPP
