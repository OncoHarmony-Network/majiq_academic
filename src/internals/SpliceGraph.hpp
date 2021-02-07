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
#include "GeneJunctions.hpp"
#include "Exons.hpp"
#include "Introns.hpp"
#include "OverGenes.hpp"


namespace majiq {
class SpliceGraph {
 protected:
  std::shared_ptr<Contigs> contigs_;
  std::shared_ptr<Genes> genes_;
  std::shared_ptr<Exons> exons_;
  std::shared_ptr<GeneJunctions> junctions_;
  std::shared_ptr<Introns> introns_;

  // indexing from contigs to overlapping genes
  std::shared_ptr<OverGenes> overgenes_;
  // indexing from genes to exon/junction/intron features
  std::vector<size_t> gene_exons_idx_;
  std::vector<size_t> gene_junctions_idx_;
  std::vector<size_t> gene_introns_idx_;

  /**
   * Index into regions that are sorted by gene
   */
  template <class ContainerT>
  static std::vector<size_t> index_gene_regions(
      const std::shared_ptr<Genes>& genes, ContainerT regions);

 public:
  // access non const pointers for use by pybind11 interface...
  std::shared_ptr<Contigs> contigs() { return contigs_; }
  std::shared_ptr<Genes> genes() { return genes_; }
  std::shared_ptr<OverGenes> overgenes() { return overgenes_; }

  std::shared_ptr<Exons> exons() { return exons_; }
  const std::vector<size_t>& gene_exons_idx() { return gene_exons_idx_; }

  std::shared_ptr<GeneJunctions> junctions() { return junctions_; }
  const std::vector<size_t>& gene_junctions_idx() {
    return gene_junctions_idx_;
  }

  std::shared_ptr<Introns> introns() { return introns_; }
  const std::vector<size_t>& gene_introns_idx() { return gene_introns_idx_; }

  // constructors
  SpliceGraph(
      const std::shared_ptr<Contigs>& contigs,
      const std::shared_ptr<Genes>& genes,
      const std::shared_ptr<Exons>& exons,
      const std::shared_ptr<GeneJunctions>& junctions,
      const std::shared_ptr<Introns>& introns)
      : contigs_{contigs},
        genes_{genes},
        exons_{exons},
        junctions_{junctions},
        introns_{introns},
        overgenes_{std::make_shared<OverGenes>(contigs_, genes_)},
        gene_exons_idx_{index_gene_regions(genes_, *exons_)},
        gene_junctions_idx_{index_gene_regions(genes_, *junctions_)},
        gene_introns_idx_{index_gene_regions(genes_, *introns_)} {
  }
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
        introns_{introns},
        overgenes_{std::make_shared<OverGenes>(contigs_, genes_)},
        gene_exons_idx_{index_gene_regions(genes_, *exons_)},
        gene_junctions_idx_{index_gene_regions(genes_, *junctions_)},
        gene_introns_idx_{index_gene_regions(genes_, *introns_)} {
  }
  SpliceGraph(const SpliceGraph& sg) = default;
  SpliceGraph(SpliceGraph&& sg) = default;
  SpliceGraph& operator=(const SpliceGraph& sg) = default;
  SpliceGraph& operator=(SpliceGraph&& sg) = default;

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


// how to index regions that are sorted by gene
template <class ContainerT>
std::vector<size_t> SpliceGraph::index_gene_regions(
    const std::shared_ptr<Genes>& genes, ContainerT regions) {
  // initialize result
  std::vector<size_t> result(genes->size() + 1);
  // index into genes
  size_t gene_idx = 0;
  for (size_t region_idx = 0; region_idx < regions.size(); ++region_idx) {
    const size_t region_gene_idx = regions[region_idx].gene.idx_;
    while (gene_idx < region_gene_idx) {
      result[++gene_idx] = region_idx;
    }
  }
  while (gene_idx < genes->size()) {
    result[++gene_idx] = regions.size();
  }
  return result;
}
}  // namespace majiq

#endif  // MAJIQ_SPLICEGRAPH_HPP
