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

// Genes
#include <algorithm>
#include <memory>
#include <vector>
#include <map>
#include <optional>
#include <boost/functional/hash.hpp>

#include "ContigRegion.hpp"
#include "Gene.hpp"
#include "MajiqTypes.hpp"



namespace majiq {
class Genes;
struct KnownGene {
 public:
  // fields
  size_t gene_idx;
  std::shared_ptr<Genes> known_genes;

  // constructors
  KnownGene() = default;
  KnownGene(size_t _gene_idx, std::shared_ptr<Genes> _known_genes)
      : gene_idx{_gene_idx}, known_genes{_known_genes} {
  }
  KnownGene(const KnownGene& x) = default;
  KnownGene(KnownGene&& x) = default;
  KnownGene& operator=(const KnownGene& x) = default;
  KnownGene& operator=(KnownGene&& x) = default;

  /**
   * Get underying gene
   */
  inline Gene& get() const;
  /**
   * Get matching known gene from another Genes object
   */
  inline KnownGene remapped(
      const std::shared_ptr<Genes>& new_known_genes) const;
  KnownContig& contig() const { return get().contig; }
  GeneStrandness& strand() const { return get().strand; }
};
// equality of known genes
inline bool operator==(const KnownGene& x, const KnownGene& y) noexcept {
  return std::tie(x.gene_idx, x.known_genes)
    == std::tie(y.gene_idx, y.known_genes);
}
// sorting/ordering based on underlying gene
inline bool operator<(const KnownGene& x, const KnownGene& y) noexcept {
  return std::tie(x.gene_idx, x.known_genes)
    < std::tie(y.gene_idx, y.known_genes);
}
// comparisons against objects with KnownGene gene or gene()
template <typename T>
inline bool operator<(const T& x, const KnownGene& y) noexcept {
  constexpr bool has_field = detail::has_gene_field<T>::value;
  constexpr bool has_function = detail::has_gene_function<T>::value;
  static_assert(has_field || has_function,
      "Type T does not have gene to compare to KnownGene");
  if (has_field) {
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
  if (has_field) {
    return x < y.gene;
  } else {
    return x < y.gene();
  }
}


class Genes : public std::enable_shared_from_this<Genes> {
 public:
  using vecGene = std::vector<Gene>;
  using mapGeneIdx = std::map<geneid_t, size_t>;

 private:
  vecGene genes_vec_;
  mapGeneIdx id_idx_map_;
  bool is_sorted_;

  static mapGeneIdx _geneid_indexes(const vecGene& v) {
    mapGeneIdx result{};
    for (size_t idx = 0; idx < v.size(); ++idx) {
      result[v[idx].gene_id()] = idx;
    }
    return result;
  }

 public:
  // access vector<string> objects
  const std::vector<geneid_t> geneids() const {
    std::vector<geneid_t> result{size()};
    std::transform(genes_vec_.begin(), genes_vec_.end(), result.begin(),
        [](const Gene& g) { return g.gene_id(); });
    return result;
  }
  const std::vector<genename_t> genenames() const {
    std::vector<genename_t> result{size()};
    std::transform(genes_vec_.begin(), genes_vec_.end(), result.begin(),
        [](const Gene& g) { return g.gene_name(); });
    return result;
  }

  Gene& get(size_t idx) { return genes_vec_[idx]; }
  KnownGene operator[](size_t idx) { return {idx, shared_from_this()}; }

  // check/add genes
  size_t size() const { return genes_vec_.size(); }
  bool contains(const geneid_t& id) const { return id_idx_map_.count(id) > 0; }
  bool contains(const Gene& x) const { return contains(x.gene_id()); }
  // get gene_idx (throw error if not present)
  size_t get_gene_idx(const geneid_t& id) const { return id_idx_map_.at(id); }
  size_t get_gene_idx(const Gene& x) const { return get_gene_idx(x.gene_id()); }
  // get gene_idx (-1 if not present)
  std::optional<size_t> safe_gene_idx(const geneid_t& id) const {
    std::optional<size_t> result;
    auto match = id_idx_map_.find(id);
    if (match != id_idx_map_.end()) {
      result = match->second;
    }
    return result;
  }
  std::optional<size_t> safe_gene_idx(const Gene& x) const {
    return safe_gene_idx(x.gene_id());
  }
  /**
   * get gene_idx (add if not present)
   */
  size_t add(const Gene& x) {
    const size_t prev_size = size();
    if (!is_sorted_ || (prev_size > 0 && genes_vec_.back() >= x)) {
      // not sorted or doesn't belong at end, so find if it is here
      auto match = id_idx_map_.find(x.gene_id());
      if (match != id_idx_map_.end()) {
        return match->second;
      } else {
        // we don't have it, and adding it will make it not sorted
        is_sorted_ = false;
      }
    }
    // not present, so we add to end
    id_idx_map_[x.gene_id()] = prev_size;
    genes_vec_.push_back(x);
    return prev_size;
  }

  /**
   * get known gene from input (combine add, operator[], requires shared_ptr
   */
  KnownGene make_known(const Gene& x) { return operator[](add(x)); }
  /**
   * get known gene from gene id (fails if doesn't exist)
   */
  KnownGene known(const geneid_t& id) { return operator[](get_gene_idx(id)); }

  Genes() : genes_vec_{}, id_idx_map_{}, is_sorted_{true} {}
  explicit Genes(const vecGene& genes)
      : genes_vec_{genes},
        id_idx_map_{_geneid_indexes(genes_vec_)},
        is_sorted_{std::is_sorted(genes_vec_.begin(), genes_vec_.end())} {
  }
  explicit Genes(vecGene&& genes)
      : genes_vec_{genes},
        id_idx_map_{_geneid_indexes(genes_vec_)},
        is_sorted_{std::is_sorted(genes_vec_.begin(), genes_vec_.end())} {
  }
  Genes(const Genes& x)
      : genes_vec_{x.genes_vec_}, id_idx_map_{}, is_sorted_{true} {
    std::sort(genes_vec_.begin(), genes_vec_.end());
    id_idx_map_ = _geneid_indexes(genes_vec_);
  }
  Genes(Genes&& x) : genes_vec_{x.genes_vec_}, id_idx_map_{}, is_sorted_{true} {
    std::sort(genes_vec_.begin(), genes_vec_.end());
    id_idx_map_ = _geneid_indexes(genes_vec_);
  }
  /**
   * Shared pointer to Genes object with same genes but guaranteed sorted order
   * (if sorted, return self, otherwise create sorted copy)
   */
  std::shared_ptr<Genes> sorted() {
    return is_sorted_ ? shared_from_this() : std::make_shared<Genes>(*this);
  }
  std::shared_ptr<const Genes> sorted() const {
    return is_sorted_ ? shared_from_this() : std::make_shared<Genes>(*this);
  }
  // view into genes object for pybind11
  const std::vector<Gene>& data() { return genes_vec_; }
  // expose whether it is sorted
  inline bool is_sorted() const noexcept { return is_sorted_; }

  // equality of genes
  friend inline bool operator==(const Genes& x, const Genes&y) noexcept {
    return x.genes_vec_ == y.genes_vec_;
  }
};

// implement KnownGene::remapped and KnownGene::get using Genes definition
inline Gene& KnownGene::get() const { return known_genes->get(gene_idx); }
inline KnownGene KnownGene::remapped(
    const std::shared_ptr<Genes>& new_known_genes) const {
  // get idx for this gene in new_known_genes
  size_t new_gene_idx = new_known_genes->add(get());
  return {new_gene_idx, new_known_genes};
}

// specialize boost::hash_value for KnownGene
inline std::size_t hash_value(const KnownGene& x) {
  std::size_t result = std::hash<size_t>{}(x.gene_idx);
  boost::hash_combine(result, x.known_genes);
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
