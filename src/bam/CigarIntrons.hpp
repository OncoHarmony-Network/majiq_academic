/**
 * CigarIntrons.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_BAM_CIGARINTRONS_HPP
#define MAJIQ_BAM_CIGARINTRONS_HPP

#include <algorithm>
#include <vector>
#include <utility>

#include "htslib/sam.h"
#include "../MajiqTypes.hpp"
#include "CigarUtils.hpp"

namespace majiq
namespace bam {
using aln_off_t = int64_t;

template<uint32_t min_overhang>
class CigarIntrons : CigarBase {
 private:
  // junctions present in the alignment
  std::vector<std::pair<position_t, position_t>> junctions_;
  // genomic/read offsets to determine offset of introns
  std::vector<std::pair<position_t, aln_off_t>> genomic_alignment_offsets_;
  // first valid value for intron end
  position_t first_pos_;
  // last valid value for intron start
  position_t last_pos_;

 public:
  /**
   * First valid value for intron end
   */
  position_t first_pos() const { return first_pos_; }
  /**
   * Last valid value for intron start
   */
  position_t last_pos() const { return last_pos_; }

  /**
   * Does the read include a junction intersecting closed interval?
   *
   * @param start, end: start and end coordinates of closed interval
   */
  bool has_junction_intron(position_t start, position_t end) const {
    for (auto&& [j_start, j_end] : junctions_) {
      if (j_start <= end && start <= j_end) {
        // if junction interval intersects intron interval
        return true;
      } else if (j_end > end) {
        // junctions_ is sorted and non-self-intersecting
        // all future junctions will come after query intron
        return false;
      }
    }
    // no junctions intersected
    return false;
  }

  /**
   * Read offset of the input position
   *
   * @param c_offset coordinate offset for open/closed query position
   * relative to alignment
   * @param query genomic coordinate to obtain offset for on aligned
   * read 
   */
  template<int64_t c_offset>
  int64_t read_offset(position_t ref_pos, bool is_start) const {
    // where is start of alignment vs reference position?
    int64_t relative_offset = (c_offset + genomic_pos_) - ref_pos;
    if (relative_offset < 0) {
      // we want offset on the alignment, not reference sequence,
      // so we need to use genomic_alignment_offsets_.
      // we want last offset before given genomic position, so
      // find first one past...
      auto first_after = std::find_if(genomic_alignment_offsets_.begin() + 1,
          genomic_alignment_offsets_.end(),
          [&ref_pos](auto x) -> bool { return x.first >= ref_pos; });
      // so last before precedes first_after
      auto last_before = first_after - 1;  // preceding iterator
      // relative offset is...
      relative_offset =
        // last genomic coordinate before ref_pos as
        // c_offset-indexed position
        (c_offset + last_before->first)
        // relative to reference position (which comes after,
        // so this is negative)
        - ref_pos  // assumed to be (c_offset-indexed)
        // adjust for where this is on the alignment
        - last_before->second;
    }
    return
      // offset relative to alignment
      relative_offset + alignment_length_ - min_overhang
      // make into offset relative to *read* --> adjust for soft
      // clipping on right
      + clipping_lengths_.second;
  }

  CigarIntrons(
      position_t genomic_pos, int32_t read_length,
      uint32_t* cigar, uint32_t n_cigar)
      : CigarBase(genomic_pos, read_length, cigar, n_cigar) {
    position_t g_pos = genomic_pos_;
    aln_off_t alignment_offset = 0;
    // initial genomic_alignment pair
    genomic_alignment_offsets_.push_back({g_pos, alignment_offset});

    // parse cigar operations to fill junctions_ and genomic_alignment_offsets_
    for (size_t i = 0; i < n_cigar_; ++i) {
      const char cigar_op = bam_cigar_op(cigar_[i]);
      const uint32_t cigar_oplen = bam_cigar_oplen(cigar_[i]);
      // if junction, add it
      if (cigar_op == BAM_CREF_SKIP) {
        junctions_.push_back({g_pos, g_pos + cigar_oplen + 1});
      }
      // advancing read or reference?
      const char cigar_type = bam_cigar_type(cigar_op);
      const bool advance_alignment
        = cigar_type & details::CIGAR_CONSUMES_QUERY;
      const bool advance_reference
        = cigar_type & details::CIGAR_CONSUMES_REFERENCE;
      // does current operation advance alignment?
      if (advance_alignment) {
        // does this operation give us the first position for valid feature?
        if ((alignment_offset < min_overhang)
            && (alignment_offset + cigar_oplen >= min_overhang)) {
          first_pos_
            = g_pos + (advance_reference ? min_overhang - alignment_offset : 0);
        }
        if ((alignment_offset + cigar_oplen > alignment_length_ - min_overhang)
            && (alignment_offset <= alignment_length_ - min_overhang)) {
          // last valid value for intron_start --> need +1
          // for 1-indexed intron start
          last_pos_
            = (1 + genomic_pos_)
            + (advance_reference
                ? alignment_length_ - min_overhang - alignment_offset
                : 0);
        }
        alignment_offset += cigar_oplen;
      }
      if (advance_reference) {
        g_pos += cigar_oplen;
      }
      if (advance_reference != advance_alignment) {
        genomic_alignment_offsets_.push_back({g_pos, alignment_offset});
      }
    }
  }
};
}  // namespace bam
}  // namespace majiq

#endif  // MAJIQ_BAM_CIGARINTRONS_HPP
