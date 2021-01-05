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
#include <boost/functional/hash.hpp>

#include "Contigs.hpp"
#include "Interval.hpp"


namespace majiq {

typedef std::string geneid_t;
typedef std::string genename_t;

enum class GeneStrandness : char {
  FORWARD = '+',
  REVERSE = '-',
  AMBIGUOUS = '.',
};
std::ostream& operator<<(std::ostream& os, const GeneStrandness& x) noexcept {
  os << static_cast<char>(x);
  return os;
}

struct Gene {
 public:
  // location
  KnownContig contig;
  ClosedInterval interval;
  GeneStrandness strand;
  // identification
  geneid_t geneid;  // unique key
  genename_t genename;

  // constructors
  Gene(KnownContig _contig, ClosedInterval _interval, GeneStrandness _strand,
      geneid_t _geneid, genename_t _genename)
      : contig{_contig},
        interval{_interval},
        strand{_strand},
        geneid{_geneid},
        genename{_genename} {
  }
  Gene(const Gene& g)
      : Gene(g.contig, g.interval, g.strand, g.geneid, g.genename) {
  }
};
// ordering with respect to genomic position
inline bool operator<(const Gene& lhs, const Gene& rhs) {
  return std::tie(lhs.contig, lhs.interval, lhs.strand, lhs.geneid)
    < std::tie(rhs.contig, rhs.interval, rhs.strand, rhs.geneid);
}
// equality only by gene id
inline bool operator==(const Gene& lhs, const Gene& rhs) {
  return lhs.geneid == rhs.geneid;
}
// allow Gene to be passed into output stream (e.g. std::cout)
std::ostream& operator<<(std::ostream& os, const Gene& x) noexcept {
  os << x.geneid;
  return os;
}

// comparison to contigs
inline bool operator<(const Gene& x, const KnownContig& y) {
  return x.contig < y;
}
inline bool operator<(const KnownContig& x, const Gene& y) {
  return x < y.contig;
}

// derived comparisons (Gene, Gene)
inline bool operator>(const Gene& x, const Gene& y) { return y < x; }
inline bool operator<=(const Gene& x, const Gene& y) { return !(y < x); }
inline bool operator>=(const Gene& x, const Gene& y) { return !(x < y); }
inline bool operator!=(const Gene& x, const Gene& y) { return !(x == y); }

// derived comparisons (Gene, KnownContig)
inline bool operator>(const KnownContig& x, const Gene& y) noexcept {
  return y < x;
}
inline bool operator>(const Gene& x, const KnownContig y) noexcept {
  return y < x;
}
inline bool operator<=(const KnownContig& x, const Gene& y) noexcept {
  return !(y < x);
}
inline bool operator<=(const Gene& x, const KnownContig y) noexcept {
  return !(y < x);
}
inline bool operator>=(const KnownContig& x, const Gene& y) noexcept {
  return !(x < y);
}
inline bool operator>=(const Gene& x, const KnownContig y) noexcept {
  return !(x < y);
}

// override boost::hash
std::size_t hash_value(const GeneStrandness& x) noexcept {
  return std::hash<char>{}(static_cast<char>(x));
}
std::size_t hash_value(const Gene& x) noexcept {
  return std::hash<geneid_t>{}(x.geneid);
}
}  // namespace majiq
// override std::hash
namespace std {
template <> struct hash<majiq::GeneStrandness> {
  std::size_t operator()(const majiq::GeneStrandness& x) const noexcept {
    return majiq::hash_value(x);
  }
};
template <> struct hash<majiq::Gene> {
  std::size_t operator()(const majiq::Gene& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
class Genes;
struct KnownGene {
 public:
  // fields
  size_t gene_idx;
  std::shared_ptr<Genes> known_genes;
  // constructors
  KnownGene(size_t _gene_idx, std::shared_ptr<Genes> _known_genes)
      : gene_idx{_gene_idx}, known_genes{_known_genes} {
  }
  KnownGene(const KnownGene& x) : KnownGene(x.gene_idx, x.known_genes) { }
  /**
   * Get underying gene
   */
  Gene& get() const;
  /**
   * Get matching known gene from another Genes object
   */
  KnownGene remapped(const std::shared_ptr<Genes>& new_known_genes) const;
  // sorting/equality based on underlying gene
  bool operator<(const KnownGene& rhs) const {
    return gene_idx < rhs.gene_idx;
  }
  bool operator==(const KnownGene& rhs) const {
    return gene_idx == rhs.gene_idx && known_genes == rhs.known_genes;
  }
};

class Genes : public std::enable_shared_from_this<Genes> {
 public:
  using vecGene = std::vector<Gene>;
  using mapGeneIdx = std::map<geneid_t, size_t>;

 private:
  vecGene genes_vec_;
  mapGeneIdx id_idx_map_;
  bool is_sorted_;

  // static functions to enable constructors to always be sorted
  template <typename GeneIt>
  static vecGene _sorted_genes_vector(GeneIt start, GeneIt end) {
    static_assert(
        std::is_same<Gene,
          typename std::iterator_traits<GeneIt>::value_type>::value,
        "GeneIt must have Gene values in _sorted_genes_vector()");
    vecGene result{start, end};
    std::sort(result.begin(), result.end());
    return result;
  }
  static mapGeneIdx _geneid_indexes(vecGene v) {
    mapGeneIdx result{};
    for (size_t idx = 0; idx < v.size(); ++idx) {
      result[v[idx].geneid] = idx;
    }
    return result;
  }

 public:
  Gene& get(size_t idx) { return genes_vec_[idx]; }
  KnownGene operator[](size_t idx) { return {idx, shared_from_this()}; }

  // check/add genes
  size_t size() const { return genes_vec_.size(); }
  bool contains(const geneid_t& id) const { return id_idx_map_.count(id) > 0; }
  bool contains(const Gene& x) const { return contains(x.geneid); }
  // get gene_idx (throw error if not present)
  size_t get_gene_idx(const geneid_t& id) const { return id_idx_map_.at(id); }
  size_t get_gene_idx(const Gene& x) const { return get_gene_idx(x.geneid); }
  /**
   * get gene_idx (add if not present)
   */
  size_t add(const Gene& x) {
    const size_t prev_size = size();
    if (!is_sorted_ || (prev_size > 0 && genes_vec_.back() >= x)) {
      // not sorted or doesn't belong at end, so find if it is here
      auto match = id_idx_map_.find(x.geneid);
      if (match != id_idx_map_.end()) {
        return match->second;
      } else {
        // we don't have it, and adding it will make it not sorted
        is_sorted_ = false;
      }
    }
    // not present, so we add to end
    id_idx_map_[x.geneid] = prev_size;
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

  // constructors always ensure that result is sorted
  Genes() : genes_vec_{}, id_idx_map_{}, is_sorted_{true} {}
  template <class GeneIt>
  Genes(GeneIt start, GeneIt end)
      : genes_vec_{_sorted_genes_vector(start, end)},
        id_idx_map_{_geneid_indexes(genes_vec_)},
        is_sorted_{true} {
  }
  Genes(const Genes& x) : Genes{x.genes_vec_.begin(), x.genes_vec_.end()} { }
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
};

// implement KnownGene::remapped and KnownGene::get using Genes definition
Gene& KnownGene::get() const { return known_genes->get(gene_idx); }
KnownGene KnownGene::remapped(
    const std::shared_ptr<Genes>& new_known_genes) const {
  // get idx for this gene in new_known_genes
  size_t new_gene_idx = new_known_genes->add(get());
  return {new_gene_idx, new_known_genes};
}

// specialize boost::hash_value for KnownGene
std::size_t hash_value(const KnownGene& x) {
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
