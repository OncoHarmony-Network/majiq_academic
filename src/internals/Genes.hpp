/**
 * Genes.hpp
 *
 * Genes for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_GENES_HPP
#define MAJIQ_GENES_HPP

#include <iostream>
#include <string>
#include <tuple>
#include <utility>

// Genes
#include <algorithm>
#include <memory>
#include <vector>
#include <map>
#include <optional>
#include <boost/functional/hash.hpp>

#include "ContigRegion.hpp"
#include "Regions.hpp"
#include "Meta.hpp"
#include "KnownFeatures.hpp"
#include "Gene.hpp"
#include "MajiqTypes.hpp"



namespace majiq {
class Genes;
using KnownGene = detail::KnownFeature<Genes>;

// comparisons against objects with KnownGene gene or gene()
template <typename T>
inline bool operator<(const T& x, const KnownGene& y) noexcept {
  constexpr bool has_field = detail::has_gene_field<T>::value;
  constexpr bool has_function = detail::has_gene_function<T>::value;
  static_assert(has_field || has_function,
      "Type T does not have gene to compare to KnownGene");
  if constexpr(has_field) {
    return x.gene < y;
  } else {
    return x.gene() < y;
  }
}
template <typename T>
inline bool operator<(const KnownGene& x, const T& y) noexcept {
  constexpr bool has_field = detail::has_gene_field<T>::value;
  constexpr bool has_function = detail::has_gene_function<T>::value;
  static_assert(has_field || has_function,
      "Type T does not have gene to compare to KnownGene");
  if constexpr(has_field) {
    return x < y.gene;
  } else {
    return x < y.gene();
  }
}


class Genes
    : public detail::KnownFeatures<detail::Regions<Gene, true>>,
      public std::enable_shared_from_this<Genes> {
 public:
  // implement creation/viewing of KnownGene as defined
  const KnownGene operator[](size_t idx) { return KnownGene{idx, shared_from_this()}; }
  // access vector<string> objects
  const std::vector<geneid_t> geneids() const {
    std::vector<geneid_t> result{size()};
    std::transform(features_.begin(), features_.end(), result.begin(),
        [](const Gene& g) { return g.gene_id(); });
    return result;
  }
  const std::vector<genename_t> genenames() const {
    std::vector<genename_t> result{size()};
    std::transform(features_.begin(), features_.end(), result.begin(),
        [](const Gene& g) { return g.gene_name(); });
    return result;
  }
  // view into genes object for pybind11
  const std::vector<Gene>& data() { return features_.elements_; }

  explicit Genes(std::vector<Gene>&& x)
      : detail::KnownFeatures<detail::Regions<Gene, true>>(std::move(x)) { }
};

// specialize boost::hash_value for KnownGene
inline std::size_t hash_value(const KnownGene& x) {
  std::size_t result = std::hash<size_t>{}(x.idx_);
  boost::hash_combine(result, x.ptr_);
  return result;
}
}  // namespace majiq
// specialize std::hash for KnownGene
namespace std {
template <> struct hash<majiq::KnownGene> {
  std::size_t operator()(const majiq::KnownGene& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

#endif  // MAJIQ_GENES_HPP
