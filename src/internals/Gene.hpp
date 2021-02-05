/**
 * Gene.hpp
 *
 * Individual gene definition
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_GENE_HPP
#define MAJIQ_GENE_HPP

#include <iostream>
#include <string>
#include <tuple>
#include <functional>
#include <boost/functional/hash.hpp>

#include "ContigRegion.hpp"
#include "Contigs.hpp"
#include "MajiqTypes.hpp"

namespace majiq {
struct GeneInfo {
  geneid_t gene_id;  // unique key
  genename_t gene_name;
};

struct Gene : public detail::ContigRegion<ClosedInterval, GeneInfo> {
  using BaseT = detail::ContigRegion<ClosedInterval, GeneInfo>;
  geneid_t& gene_id() noexcept { return data.gene_id; }
  const geneid_t& gene_id() const noexcept { return data.gene_id; }
  genename_t& gene_name() noexcept { return data.gene_name; }
  const genename_t& gene_name() const noexcept { return data.gene_name; }
  geneid_t unique_key() const noexcept { return gene_id(); }  // for KnownFeatures
  Gene(KnownContig _contig, ClosedInterval _coordinates, GeneStrandness _strand,
      GeneInfo _info)
      : BaseT{_contig, _coordinates, _strand, _info} { }
  Gene(KnownContig _contig, ClosedInterval _coordinates, GeneStrandness _strand,
      geneid_t _geneid, genename_t _genename)
      : Gene{_contig, _coordinates, _strand, GeneInfo{_geneid, _genename}} { }
  Gene() = default;
  Gene(const Gene&) = default;
  Gene(Gene&&) = default;
  Gene& operator=(const Gene&) = default;
  Gene& operator=(Gene&&) = default;
};
// ordering with respect to genomic position
inline bool operator<(const Gene& lhs, const Gene& rhs) noexcept {
  return std::tie(lhs.contig, lhs.coordinates, lhs.strand, lhs.gene_id())
    < std::tie(rhs.contig, rhs.coordinates, rhs.strand, rhs.gene_id());
}
// equality only by gene id
inline bool operator==(const Gene& lhs, const Gene& rhs) noexcept {
  return lhs.gene_id() == rhs.gene_id();
}
// allow Gene to be passed into output stream (e.g. std::cout)
inline std::ostream& operator<<(std::ostream& os, const Gene& x) noexcept {
  os << x.gene_id();
  return os;
}

// comparison to contigs
inline bool operator<(const Gene& x, const KnownContig& y) noexcept {
  return x.contig < y;
}
inline bool operator<(const KnownContig& x, const Gene& y) noexcept {
  return x < y.contig;
}

// derived comparisons (Gene, Gene)
inline bool operator>(const Gene& x, const Gene& y) noexcept { return y < x; }
inline bool operator<=(const Gene& x, const Gene& y) noexcept {
  return !(y < x);
}
inline bool operator>=(const Gene& x, const Gene& y) noexcept {
  return !(x < y);
}
inline bool operator!=(const Gene& x, const Gene& y) noexcept {
  return !(x == y);
}

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
inline std::size_t hash_value(const Gene& x) noexcept {
  return std::hash<geneid_t>{}(x.gene_id());
}
}  // namespace majiq

namespace std {
template <> struct hash<majiq::Gene> {
  std::size_t operator()(const majiq::Gene& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

#endif  // MAJIQ_GENE_HPP
