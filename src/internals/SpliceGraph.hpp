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
#include "GeneIntrons.hpp"
#include "PassedJunctions.hpp"
#include "ExonConnections.hpp"


namespace majiq {
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

  static Exons InferExons(const Exons& source, const GeneJunctions& junctions);

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
}  // namespace majiq

#endif  // MAJIQ_SPLICEGRAPH_HPP
