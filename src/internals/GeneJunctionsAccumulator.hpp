/**
 * GeneJunctionsAccumulator.hpp
 *
 * Junctions for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_GENEJUNCTIONSACCUMULATOR_HPP
#define MAJIQ_GENEJUNCTIONSACCUMULATOR_HPP

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "Genes.hpp"
#include "GeneJunctions.hpp"

namespace majiq {

class GeneJunctionsAccumulator {

 private:
  std::shared_ptr<Genes> genes_;
  std::vector<std::map<OpenInterval, detail::ConnectionData>> gene_junctions_;

 public:
  // constructors should set genes, set gene_junctions to correct size
  explicit GeneJunctionsAccumulator(const std::shared_ptr<Genes>& genes)
      : genes_{genes}, gene_junctions_(genes->size()) {
    if (genes_ == nullptr) {
      throw std::runtime_error("GeneJunctionsAccumulator given null genes");
    }
  }
  GeneJunctionsAccumulator() = delete;
  GeneJunctionsAccumulator(const GeneJunctionsAccumulator&) = default;
  GeneJunctionsAccumulator(GeneJunctionsAccumulator&&) = default;
  GeneJunctionsAccumulator& operator=(const GeneJunctionsAccumulator&) = default;
  GeneJunctionsAccumulator& operator=(GeneJunctionsAccumulator&&) = default;

  const std::shared_ptr<Genes> genes() const { return genes_; }

  // add GeneJunctions as denovo or as annotated
  void Add(const GeneJunctions& junctions, const bool make_annotated) {
    // check that genes match
    if (genes_ != junctions.parents()) {
      throw std::invalid_argument(
          "junction genes do not match in GeneJunctionsAccumulator::Add()");
    }
    // for each junction, add or update associated value in gene_junctions_
    for (auto&& x : junctions) {
      // map from intervals to junction flags for gene associated with junction
      auto& junction_map = gene_junctions_[x.gene.idx_];
      // get data from x, updating denovo to false if make_annotated set
      detail::ConnectionData x_data = x.data;
      if (make_annotated) { x_data.denovo = false; }
      // try to emplace new junction into map.
      // If it exists, get iterator to it (but don't update it yet)
      auto [pair_it, new_junction] =
          junction_map.try_emplace(x.coordinates, x_data);
      // If the junction existed already, update flags
      if (!new_junction) {
        // get reference to existing junction data
        detail::ConnectionData& map_data = pair_it->second;
        // update flags: denovo or simplified only if both are, passed if either
        map_data.denovo = x_data.denovo && map_data.denovo;
        map_data.passed_build = x_data.passed_build || map_data.passed_build;
        map_data.simplified = x_data.simplified && map_data.simplified;
      }
    }
    return;
  }

  // return GeneJunctions with accumulated junctions
  GeneJunctions Accumulated() const {
    std::vector<GeneJunction> junctions_vec;
    // add junctions per gene in sorted order
    for (size_t gene_idx = 0; gene_idx < gene_junctions_.size(); ++gene_idx) {
      for (const auto& [coordinates, data] : gene_junctions_[gene_idx]) {
        junctions_vec.emplace_back(KnownGene(gene_idx, genes_), coordinates,
            data.denovo, data.passed_build, data.simplified);
      }
    }
    return GeneJunctions(genes_, std::move(junctions_vec));
  }
};

}  // namespace majiq


#endif  // MAJIQ_GENEJUNCTIONSACCUMULATOR_HPP
