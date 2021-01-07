/**
 * Contigs.hpp
 *
 * Contigs (i.e. chromosomes) for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_CONTIGS_HPP
#define MAJIQ_CONTIGS_HPP

#include <iterator>
#include <functional>
#include <string>
#include <iostream>

// Contigs
#include <memory>
#include <map>
#include <vector>
#include <algorithm>
#include <boost/functional/hash.hpp>


// define individual contig
namespace majiq {

using seqid_t = std::string;  // type for contig ids

struct Contig {
 public:
  seqid_t seqid;

  explicit Contig(const seqid_t& _seqid) : seqid{_seqid} {}
  explicit Contig(seqid_t&& _seqid) : seqid{_seqid} {}
  Contig(const Contig& x) = default;
  Contig(Contig&& x) = default;
  Contig& operator=(const Contig& x) = default;
  Contig& operator=(Contig&& x) = default;
};
inline bool operator<(const Contig& x, const Contig& y) noexcept {
  return x.seqid < y.seqid;
}
inline bool operator==(const Contig& x, const Contig& y) noexcept {
  return x.seqid == y.seqid;
}
inline std::ostream& operator<<(std::ostream& os, const Contig& x) noexcept {
  os << x.seqid;
  return os;
}
// override boost::hash<Contig>
std::size_t hash_value(const Contig& x) noexcept {
  return std::hash<seqid_t>{}(x.seqid);
}
}  // namespace majiq

namespace std {
// override std::hash<Contig>
template <> struct hash<majiq::Contig> {
  std::size_t operator()(const majiq::Contig& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std

namespace majiq {
class Contigs;  // forward declaration implemented later in file

struct KnownContig {
 public:
  // data
  size_t contig_idx;
  std::shared_ptr<Contigs> known_contigs;

  // access true contig data
  const Contig& get() const;

  // constructors
  KnownContig(
      size_t _contig_idx,
      const std::shared_ptr<Contigs>& _known_contigs)
      : contig_idx{_contig_idx},
        known_contigs{_known_contigs} {
  }
  KnownContig(const KnownContig& x) = default;
  KnownContig(KnownContig&& x) = default;
  KnownContig& operator=(const KnownContig& x) = default;
  KnownContig& operator=(KnownContig&& x) = default;
};
inline bool operator<(const KnownContig& x, const KnownContig& y) noexcept {
  return x.contig_idx < y.contig_idx;
}
inline bool operator==(const KnownContig& x, const KnownContig& y) noexcept {
  return x.contig_idx == y.contig_idx;
}

class Contigs : public std::enable_shared_from_this<Contigs> {
 private:
  std::map<Contig, size_t> contig_idx_map_;
  std::vector<Contig> contigs_vec_;

 public:
  // retrieve contigs
  const Contig& get(size_t idx) const { return contigs_vec_[idx]; }
  /**
   * Access KnownContig, which requires Contigs managed by shared_ptr
   */
  const KnownContig operator[](size_t idx) {
    return KnownContig{idx, shared_from_this()};
  }

  // check/add contigs
  size_t size() const { return contigs_vec_.size(); }
  bool contains(const Contig& x) const { return contig_idx_map_.count(x) > 0; }
  bool contains(const seqid_t& x) const { return contains(Contig{x}); }
  size_t get_contig_idx(const Contig& x) const { return contig_idx_map_.at(x); }
  size_t get_contig_idx(const seqid_t& x) const {
    return get_contig_idx(Contig{x});
  }
  /**
   * Existing or new contig_idx for input contig
   */
  size_t add(const Contig& x) {
    // get iterator to match (or end if match doesn't exist)
    auto match = contig_idx_map_.find(x);
    if (match == contig_idx_map_.end()) {
      const size_t new_contig_idx = size();
      contig_idx_map_[x] = new_contig_idx;
      contigs_vec_.push_back(x);
      return new_contig_idx;
    } else {
      return match->second;
    }
  }
  size_t add(const seqid_t& x) { return add(Contig{x}); }

  /**
   * get known contig from input (combine add and operator[]
   *
   * @note: requires this to be managed by shared_ptr
   */
  const KnownContig make_known(const Contig& x) { return operator[](add(x)); }
  const KnownContig make_known(const seqid_t& x) { return operator[](add(x)); }

  // constructors
  Contigs() = default;
  Contigs(const Contigs& x) = default;
  Contigs(Contigs&& x) = default;
  Contigs& operator=(const Contigs& x) = default;
  Contigs& operator=(Contigs&& x) = default;
  template <class It>
  Contigs(It first, It last) : Contigs{} {
    using TypeIt = typename std::iterator_traits<It>::value_type;
    static_assert(
        std::is_same<Contig, TypeIt>::value
        || std::is_same<seqid_t, TypeIt>::value,
        "Iterator constructor must be given values of Contig or seqid_t");
    // add each contig from the iterator
    for (It it = first; it != last; ++it) {
      add(*it);
    }
  }

  // get vector of seqids
  const std::vector<seqid_t> seqids() const {
    std::vector<seqid_t> result{size()};
    std::transform(contigs_vec_.begin(), contigs_vec_.end(), result.begin(),
        [](const Contig& x) -> seqid_t { return x.seqid; });
    return result;
  }

  // make printable
  friend std::ostream& operator<<(std::ostream& os, const Contigs& x) noexcept;
};
std::ostream& operator<<(std::ostream& os, const Contigs& contigs) noexcept {
  os << "Contigs[";
  if (contigs.size() > 0) {
    for (size_t i = 0; i < contigs.size(); ++i) {
      os << contigs.get(i) << (i < contigs.size() - 1 ? ", " : "]");
    }
  } else {
    os << "]";
  }
  return os;
}
// implement forward declaration of get() now that we've defined Contigs
const Contig& KnownContig::get() const {
  return known_contigs->get(contig_idx);
}
// specialize hash for KnownContig
std::size_t hash_value(const KnownContig& x) noexcept {
  std::size_t result = std::hash<size_t>{}(x.contig_idx);
  boost::hash_combine(result, x.known_contigs);
  return result;
}
}  // namespace majiq
namespace std {
// override std::hash<KnownContig>
template <> struct hash<majiq::KnownContig> {
  std::size_t operator()(const majiq::KnownContig& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}

#endif  // MAJIQ_CONTIGS_HPP
