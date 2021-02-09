/**
 * Contigs.hpp
 *
 * Contigs (i.e. chromosomes) for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_CONTIGS_HPP
#define MAJIQ_CONTIGS_HPP

#include <iostream>
#include <vector>
#include <algorithm>

#include "Contig.hpp"
#include "KnownFeatures.hpp"
#include "Meta.hpp"

namespace majiq {
class Contigs;  // forward declaration to define KnownContig

using KnownContig = detail::KnownFeature<Contigs>;
// enable comparisons against objects with KnownContig contig or contig()
template <typename T,
         std::enable_if_t<detail::has_contig_field<T>::value, bool> = true>
inline bool operator<(const T& x, const KnownContig& y) noexcept {
  return x.contig < y;
}
template <typename T,
         std::enable_if_t<detail::has_contig_field<T>::value, bool> = true>
inline bool operator<(const KnownContig& x, const T& y) noexcept {
  return x < y.contig;
}
template <typename T,
         std::enable_if_t<detail::has_contig_function<T>::value, bool> = true>
inline bool operator<(const T& x, const KnownContig& y) noexcept {
  return x.contig() < y;
}
template <typename T,
         std::enable_if_t<detail::has_contig_function<T>::value, bool> = true>
inline bool operator<(const KnownContig& x, const T& y) noexcept {
  return x < y.contig();
}

class Contigs
    : public detail::KnownFeatures<std::vector<Contig>>,
      public std::enable_shared_from_this<Contigs> {
 public:
  // implement creation/viewing of KnownContig as defined
  KnownContig operator[](size_t idx) {
    return KnownContig{idx, shared_from_this()};
  }
  KnownContig begin() { return operator[](0); }
  KnownContig end() { return operator[](size()); }

  size_t add(const Contig& x) {
    auto key = x.unique_key();
    auto match = idx_map_.find(key);
    if (match == idx_map_.end()) {
      const size_t new_idx = size();
      idx_map_[key] = new_idx;
      features_.push_back(x);
      return new_idx;
    } else {
      return match->second;
    }
  }
  const KnownContig make_known(const Contig& x) { return operator[](add(x)); }

  const std::vector<seqid_t> seqids() const {
    std::vector<seqid_t> result(size());
    std::transform(features_.begin(), features_.end(), result.begin(),
        [](const Contig& x) -> seqid_t { return x.seqid; });
    return result;
  }
};
inline std::ostream& operator<<(std::ostream& os, const Contigs& x) noexcept {
  os << "Contigs[";
  if (x.size() > 0) {
    for (size_t i = 0; i < x.size(); ++i) {
      os << x.get(i) << (i < x.size() - 1 ? ", " : "]");
    }
  } else {
    os << "]";
  }
  return os;
}
}  // namespace majiq

#endif  // MAJIQ_CONTIGS_HPP
