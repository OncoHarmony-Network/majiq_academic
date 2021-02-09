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
struct OverGene;

class Genes
    : public detail::KnownFeatures<detail::Regions<Gene, true>>,
      public std::enable_shared_from_this<Genes> {
 public:
  // forward declaration of creation/viewing of KnownGene
  KnownGene operator[](size_t idx);
  KnownGene begin();
  KnownGene end();
  OverGene overgene_begin();
  OverGene overgene_end();
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
  const std::vector<position_t>& position_cummax() {
    return features_.elements_end_cummax_;
  }

  explicit Genes(std::vector<Gene>&& x)
      : detail::KnownFeatures<detail::Regions<Gene, true>>(std::move(x)) { }
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
        || coordinates().start > (*this - 1).coordinates().end);
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

// define overgenes
struct OverGeneData {
  KnownGene first_;
  KnownGene last_;

  const std::shared_ptr<Genes>& genes_ptr() const { return first_.ptr_; }
  bool valid() const noexcept {
    return genes_ptr() != nullptr && first_.idx_ < genes_ptr()->size();
  }
  const KnownContig contig() const {
    return valid() ? KnownContig{} : first_.contig();
  }
  const ClosedInterval coordinates() const {
    return valid()
      ? ClosedInterval{
          first_.coordinates().start,
          (last_ - 1).position_cummax()}
      : ClosedInterval{};
  }

  // get last known gene (default constructed if no pointer)
  explicit OverGeneData(std::shared_ptr<Genes> genes_ptr)
      : first_{genes_ptr == nullptr ? KnownGene{} : genes_ptr->end()},
        last_{genes_ptr == nullptr ? KnownGene{} : genes_ptr->end()} { }
  OverGeneData() : OverGeneData{std::shared_ptr<Genes>{}} { }
  // overgene that contains given gene
  explicit OverGeneData(KnownGene gene)
      : first_{gene.IsOverGeneStart() ? gene : gene.PrevOverGeneStart()},
        last_{gene.NextOverGeneStart()} { }

  OverGeneData& operator++() {
    std::swap(first_, last_);  // first is previous last
    last_ = first_.NextOverGeneStart();  // last updated to next overgene start
    return *this;
  }
  OverGeneData& operator--() {
    std::swap(first_, last_);  // last is previous first
    first_ = last_.PrevOverGeneStart();
    return *this;
  }
};
class OverGene : public detail::ContigRegion<ClosedInterval, OverGeneData> {
 private:
  void UpdateLocation() {
    contig = data.contig();
    coordinates = data.coordinates();
    return;
  }

 public:
  using BaseT = detail::ContigRegion<ClosedInterval, OverGeneData>;
  const KnownGene& first() const { return data.first_; }
  const KnownGene& last() const { return data.last_; }
  OverGene& operator++() { ++data; UpdateLocation(); return *this; }
  OverGene& operator--() { --data; UpdateLocation(); return *this; }

  /**
   * Overgene containing input gene
   */
  explicit OverGene(OverGeneData x)
      : BaseT{x.contig(), x.coordinates(), GeneStrandness::AMBIGUOUS, x} { }
  explicit OverGene(const KnownGene& x)
      : OverGene{OverGeneData{x}} { }
  OverGene(std::shared_ptr<Genes> genes_ptr, bool make_begin)
      : OverGene{
          make_begin
          ? OverGeneData{KnownGene{0, genes_ptr}}
          : OverGeneData{genes_ptr}} { }

  friend inline bool operator==(const OverGene& x, const OverGene& y) noexcept {
    return x.first().idx_ == y.first().idx_;
  }
  friend inline bool operator!=(const OverGene& x, const OverGene& y) noexcept {
    return !(x == y);
  }
};

inline OverGene Genes::overgene_begin() {
  return OverGene{shared_from_this(), true};
}
inline OverGene Genes::overgene_end() {
  return OverGene{shared_from_this(), false};
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
