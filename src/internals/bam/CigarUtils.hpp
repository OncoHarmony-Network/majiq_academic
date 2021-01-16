/**
 * CigarUtils.hpp
 *
 * Helper functions and base class for working with CIGAR strings
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_BAM_CIGARUTILS_HPP
#define MAJIQ_BAM_CIGARUTILS_HPP

#include <htslib/sam.h>
#include <utility>
#include "../MajiqTypes.hpp"

namespace majiq {
namespace bam {
namespace detail {
// cigar type flags
constexpr char CIGAR_CONSUMES_QUERY = 1;
constexpr char CIGAR_CONSUMES_REFERENCE = 2;

/**
 * Two tasks:
 * (1) adjust (by reference) cigar pointer, n_cigar to ignore clipping
 * (2) return length of query sequence removed on left/right due to clipping
 *
 * @param cigar array of cigar operations passed by mutable reference
 * @param n_cigar number of cigar operations passed by mutable reference
 */
std::pair<int32_t, int32_t> adjust_cigar_soft_clipping(
    uint32_t*& cigar, uint32_t& n_cigar);
}  // namespace detail

/**
 * Base class for working with CIGAR strings
 * Acquires relevant information after removing contribution from clipping
 */
class CigarBase {
 protected:
  // mapping position (genome/target) at start of alignment
  const position_t genomic_pos_;
  // read length ~ left-clipping + alignment + right-clipping
  const std::pair<int32_t, int32_t> clipping_lengths_;
  const int32_t alignment_length_;
  // cigar operation for alignment (not read, after removing clipping)
  const uint32_t* cigar_;
  const uint32_t n_cigar_;

 public:
  CigarBase(
      position_t genomic_pos,
      int32_t read_length,
      uint32_t* cigar,
      uint32_t n_cigar)
      : genomic_pos_{genomic_pos},
        clipping_lengths_{detail::adjust_cigar_soft_clipping(cigar, n_cigar)},
        alignment_length_{
          read_length - clipping_lengths_.first - clipping_lengths_.second},
        cigar_{cigar},
        n_cigar_{n_cigar} {
  }
};

// implementation of detail::adjust_cigar_soft_clipping
std::pair<int32_t, int32_t> detail::adjust_cigar_soft_clipping(
    uint32_t*& cigar, uint32_t& n_cigar) {
  // initialize lengths of soft clipping on left/right to return
  int32_t left_length = 0;
  int32_t right_length = 0;
  // remove clipping on the right
  for (uint32_t i = n_cigar - 1; i >= 0; --i) {  // iterate backwards
    const char cigar_op = bam_cigar_op(cigar[i]);
    if (cigar_op == BAM_CHARD_CLIP) {
      // ignore hard clipping cigar operations on right
      --n_cigar;
    } else if (cigar_op == BAM_CSOFT_CLIP) {
      // ignore soft clipping cigar operations, update right length
      --n_cigar;
      right_length += bam_cigar_oplen(cigar[i]);
    } else {
      break;
    }
  }
  // remove clipping on the left
  size_t lhs_clipping = 0;  // offset to apply to *cigar_ptr
  for (uint32_t i = 0; i < n_cigar; ++i) {
    const char cigar_op = bam_cigar_op(cigar[i]);
    if (cigar_op == BAM_CHARD_CLIP) {
      // ignore hard clipping cigar operations on left
      ++lhs_clipping;
    } else if (cigar_op == BAM_CSOFT_CLIP) {
      // ignore soft clipping cigar operations, update left length
      ++lhs_clipping;
      left_length += bam_cigar_oplen(cigar[i]);
    } else {
      break;
    }
  }
  if (lhs_clipping > 0) {
    n_cigar -= lhs_clipping;
    cigar = cigar + lhs_clipping;
  }
  return std::make_pair(left_length, right_length);
}
}  // namespace bam
}  // namespace majiq

#endif  // MAJIQ_BAM_CIGARUTILS_HPP
