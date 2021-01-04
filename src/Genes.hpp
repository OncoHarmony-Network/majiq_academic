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

  /**
   * Ordering with respect to genomic position
   */
  bool operator<(const Gene& rhs) const {
    return std::tie(contig, interval, strand, geneid)
      < std::tie(rhs.contig, rhs.interval, rhs.strand, rhs.geneid);
  }
  /**
   * Equality only by geneid
   */
  bool operator==(const Gene& rhs) const {
    return geneid == rhs.geneid;
  }
};
// allow Gene to be passed into output stream (e.g. std::cout)
std::ostream& operator<<(std::ostream& os, const Gene& x) noexcept {
  os << x.geneid;
  return os;
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
  KnownGene remapped(const std::shared_ptr<Genes> new_known_genes) const;
  // sorting/equality based on underlying gene
  bool operator<(const KnownGene& rhs) const {
    return get() < rhs.get();
  }
  bool operator==(const KnownGene& rhs) const {
    return get() == rhs.get();
  }
};

class Genes : std::enable_shared_from_this<Genes> {
 private:
  std::map<geneid_t, size_t> id_idx_map_;
  std::vector<Gene> genes_vec_;
  bool is_sorted_;

 public:
  Gene& get(size_t idx) { return genes_vec_[idx]; }
  KnownGene operator[](size_t idx) { return {idx, shared_from_this()}; }

  // check/add genes
  size_t size() const { return genes_vec_.size(); }
  bool contains(const geneid_t& id) const { return id_idx_map_.count(id) > 0; }
  bool contains(const Gene& x) const { return contains(x.geneid); }
  size_t get_gene_idx(const geneid_t& id) const { return id_idx_map_.at(id); }
  size_t get_gene_idx(const Gene& x) const { return get_gene_idx(x.geneid); }
  /**
   * Index (old or new) of specified gene
   */
  size_t add(const Gene& x) {
    if (contains(x)) {
      return get_gene_idx(x);
    } else {
      const size_t new_gene_idx = size();
      id_idx_map_[x.geneid] = new_gene_idx;
      is_sorted_ = is_sorted_ && genes_vec_[new_gene_idx - 1] < x;
      genes_vec_.push_back(x);
      return new_gene_idx;
    }
  }

  // constructors
  Genes() : id_idx_map_{}, genes_vec_{}, is_sorted_{true} {}
  Genes(const Genes& x)
      : id_idx_map_{x.id_idx_map_.begin(), x.id_idx_map_.end()},
        genes_vec_{x.genes_vec_.begin(), x.genes_vec_.end()},
        is_sorted_{x.is_sorted_} {
  }
  template <class It>
  Genes(It first, It last) : Genes{} {
    for (It it = first; it != last; ++it) {
      add(*it);
    }
  }
  /**
   * Shared pointer to Genes object with same genes but guaranteed sorted order
   * (if sorted, return self, otherwise create sorted copy)
   */
  std::shared_ptr<Genes> sorted() {
    if (is_sorted_) {
      return shared_from_this();
    } else {
      // create empty genes
      std::shared_ptr<Genes> result = std::make_shared<Genes>();
      // copy of genes vector to new empty genes, then sort
      result->genes_vec_ = genes_vec_;
      std::sort(result->genes_vec_.begin(), result->genes_vec_.end());
      // update result->id_idx_map_ with newly sorted genes
      for (size_t idx = 0; idx < result->genes_vec_.size(); ++idx) {
        result->id_idx_map_[result->genes_vec_[idx].geneid] = idx;
      }
      // it's sorted...
      result->is_sorted_ = true;
      return result;
    }
  }
};
// implement KnownGene::remapped and KnownGene::get using Genes definition
Gene& KnownGene::get() const { return known_genes->get(gene_idx); }
KnownGene KnownGene::remapped(
    const std::shared_ptr<Genes> new_known_genes) const {
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
