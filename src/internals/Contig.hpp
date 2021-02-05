/**
 * Contig.hpp
 *
 * Contig (i.e. chromosome) for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_CONTIG_HPP
#define MAJIQ_CONTIG_HPP

#include <iostream>
#include <functional>
#include <utility>

#include "MajiqTypes.hpp"

namespace majiq {
struct Contig {
 public:
  seqid_t seqid;

  // for KnownFeatures
  seqid_t unique_key() const noexcept { return seqid; }

  explicit Contig(const seqid_t& _seqid) : seqid{_seqid} {}
  explicit Contig(seqid_t&& _seqid) : seqid{std::move(_seqid)} {}
  Contig(const Contig&) = default;
  Contig(Contig&&) = default;
  Contig() = default;
  Contig& operator=(const Contig&) = default;
  Contig& operator=(Contig&&) = default;
};
inline bool operator<(const Contig& x, const Contig& y) noexcept {
  return x.seqid < y.seqid;
}
inline std::ostream& operator<<(std::ostream& os, const Contig& x) noexcept {
  os << x.seqid;
  return os;
}
inline std::size_t hash_value(const Contig& x) noexcept {
  return std::hash<seqid_t>{}(x.seqid);
}
}  // namespace majiq

namespace std {
template <> struct hash<majiq::Contig> {
  std::size_t operator()(const majiq::Contig& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

#endif
