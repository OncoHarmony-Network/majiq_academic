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
class KnownGene;

class Genes
    : public detail::KnownFeatures<detail::Regions<Gene, true>>,
      public std::enable_shared_from_this<Genes> {
 public:
  // forward declaration of creation/viewing of KnownGene
  KnownGene operator[](size_t idx);
  KnownGene begin();
  KnownGene end();
  KnownGene begin_contig(size_t idx);
  KnownGene end_contig(size_t idx);
  KnownGene begin_contig(const KnownContig& contig);
  KnownGene end_contig(const KnownContig& contig);
  // access contigs that are parent to this
  std::shared_ptr<Contigs> contigs() const { return features_.parents(); }
  std::shared_ptr<Contigs> parents() const { return contigs(); }
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
  const std::vector<Gene>& data() { return features_.data(); }
  const std::vector<size_t>& parent_idx_offsets() {
    return features_.parent_idx_offsets();
  }
  const std::vector<position_t>& position_cummax() {
    return features_.elements_end_cummax_;
  }

// only allow Genes to be created as shared_ptr by Genes::create()
 private:
  Genes(const std::shared_ptr<Contigs>& contigs, std::vector<Gene>&& x)
      : detail::KnownFeatures<detail::Regions<Gene, true>>(
          detail::Regions<Gene, true>{contigs, std::move(x)}) {
    if (parents() == nullptr) {
      throw std::invalid_argument("Genes requires non-null contigs");
    }
  }
  struct CreateKey {
    const std::shared_ptr<Contigs>& contigs_;
    std::vector<Gene>& x_;
    CreateKey(const std::shared_ptr<Contigs>& contigs, std::vector<Gene>& x)
        : contigs_{contigs}, x_{x} { }
  };

 public:
  // public constructor requiring private type
  explicit Genes(CreateKey key) : Genes{key.contigs_, std::move(key.x_)} { }
  // allowing std::make_shared to be used in private manner
  static std::shared_ptr<Genes> create(
      const std::shared_ptr<Contigs>& contigs, std::vector<Gene>&& x) {
    return std::make_shared<Genes>(CreateKey{contigs, x});
  }
};

class KnownGene : public detail::KnownFeature<Genes, KnownGene> {
 public:
  const KnownContig& contig() const { return get().contig; }
  const GeneStrandness& strand() const { return get().strand; }
  const ClosedInterval& coordinates() const { return get().coordinates; }
  const position_t& position_cummax() const {
    return ptr_->position_cummax()[idx_];
  }

  KnownGene(size_t idx, std::shared_ptr<Genes> ptr)
      : detail::KnownFeature<Genes, KnownGene>{idx, ptr} { }
  KnownGene() = default;
  KnownGene(const KnownGene&) = default;
  KnownGene(KnownGene&&) = default;
  KnownGene& operator=(const KnownGene&) = default;
  KnownGene& operator=(KnownGene&&) = default;

  /**
   * True when this gene does not overlap with previous gene
   */
  bool IsOverGeneStart() const {
    return (
        // null gene pointer
        ptr_ == nullptr
        // start
        || idx_ == 0
        // past the end
        || idx_ >= ptr_->size()
        // new contig
        || contig() != (*this - 1).contig()
        // coordinates past previous max
        || coordinates().start > (*this - 1).position_cummax());
  }
  KnownGene NextOverGeneStart() const {
    if (ptr_ == nullptr) { return *this; }
    if (*this >= ptr_->end()) {
      return ptr_->end();
    } else {
      KnownGene result = *this + 1;
      for (; !result.IsOverGeneStart(); ++result) { }
      return result;
    }
  }
  KnownGene PrevOverGeneStart() const {
    if (ptr_ == nullptr) { return *this; }
    if (*this <= ptr_->begin()) {
      return ptr_->begin();
    } else if (*this > ptr_->end()) {
      return ptr_->end();
    } else {
      KnownGene result = *this - 1;
      for (; !result.IsOverGeneStart(); --result) { }
      return result;
    }
  }
};

inline KnownGene Genes::operator[](size_t idx) {
  return KnownGene{idx, shared_from_this()};
}
inline KnownGene Genes::begin() { return operator[](0); }
inline KnownGene Genes::end() { return operator[](size()); }
inline KnownGene Genes::begin_contig(size_t idx) {
  return operator[](features_.parent_idx_offsets_[idx]);
}
inline KnownGene Genes::end_contig(size_t idx) { return begin_contig(1 + idx); }
inline KnownGene Genes::begin_contig(const KnownContig& contig) {
  return begin_contig(contig.idx_);
}
inline KnownGene Genes::end_contig(const KnownContig& contig) {
  return end_contig(contig.idx_);
}


// comparisons against objects with KnownGene gene or gene()
template <typename T,
         std::enable_if_t<detail::has_gene_field<T>::value, bool> = true>
inline bool operator<(const T& x, const KnownGene& y) noexcept {
  return x.gene < y;
}
template <typename T,
         std::enable_if_t<detail::has_gene_field<T>::value, bool> = true>
inline bool operator<(const KnownGene& x, const T& y) noexcept {
  return x < y.gene;
}
template <typename T,
         std::enable_if_t<detail::has_gene_function<T>::value, bool> = true>
inline bool operator<(const T& x, const KnownGene& y) noexcept {
  return x.gene() < y;
}
template <typename T,
         std::enable_if_t<detail::has_gene_function<T>::value, bool> = true>
inline bool operator<(const KnownGene& x, const T& y) noexcept {
  return x < y.gene();
}

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
