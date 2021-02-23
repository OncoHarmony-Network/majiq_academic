/**
 * MaijqTypes.hpp
 *
 * base types for MAJIQ
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_TYPES_HPP
#define MAJIQ_TYPES_HPP

#include <cstdint>
#include <string>
#include <iostream>
#include <functional>
#include <boost/functional/hash.hpp>

namespace majiq {

struct EmptyDataT {};

using real_t = float;  // for numerics
using position_t = int64_t;  // for coordinates, difference in coordinates
using junction_pos_t = int32_t;  // position for junctions
using junction_ct_t = int32_t;  // counts of reads with splits

using seqid_t = std::string;  // type for contig ids
using geneid_t = std::string;
using genename_t = std::string;

enum class GeneStrandness : unsigned char {
  FORWARD = '+',
  REVERSE = '-',
  AMBIGUOUS = '.',
};
inline std::ostream& operator<<(
    std::ostream& os, const GeneStrandness& x) noexcept {
  os << static_cast<char>(x);
  return os;
}
inline std::size_t hash_value(const GeneStrandness& x) noexcept {
  return std::hash<char>{}(static_cast<char>(x));
}

enum class ExperimentStrandness : unsigned char {
  FORWARD = 'F',  // read 1 is in forward direction, Salmon ISF
  REVERSE = 'R',  // read 2 is in forward direction, Salmon ISR
  NONE = 'N',  // could be either way, Salmon IU
};
}  // namespace majiq

namespace std {
template <> struct hash<majiq::GeneStrandness> {
  std::size_t operator()(const majiq::GeneStrandness& x) const noexcept {
    return majiq::hash_value(x);
  }
};
}  // namespace std


#endif  // MAJIQ_TYPES_HPP
