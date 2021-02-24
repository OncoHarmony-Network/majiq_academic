/**
 * CigarRegions.hpp
 *
 * Iterate over aligned/unaligned regions for a read working with CIGAR strings
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_BAM_CIGARREGIONS_HPP
#define MAJIQ_BAM_CIGARREGIONS_HPP

#include <htslib/sam.h>

#include <algorithm>
#include <utility>

#include "../MajiqTypes.hpp"
#include "../Interval.hpp"

namespace majiq {
namespace bam {

class CigarRegions {
 public:
  // mapping position at start of alignment
  const position_t genomic_pos_;
  // read length ~ left-clipping + alignment + right-clipping
  const std::pair<int32_t, int32_t> clipping_lengths_;
  const int32_t alignment_length_;

 private:
  // cigar operations
  const uint32_t* cigar_;
  const uint32_t n_cigar_;

  // static const values/helper functions
  static constexpr char CIGAR_CONSUMES_QUERY = 1;
  static constexpr char CIGAR_CONSUMES_REFERENCE = 2;
  static std::pair<int32_t, int32_t> adjust_cigar_soft_clipping(
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

 public:
  CigarRegions(
      position_t genomic_pos, int32_t read_length,
      uint32_t* cigar, uint32_t n_cigar)
      : genomic_pos_{genomic_pos},
        clipping_lengths_{adjust_cigar_soft_clipping(cigar, n_cigar)},
        alignment_length_{
          read_length - clipping_lengths_.first - clipping_lengths_.second},
        cigar_{cigar},
        n_cigar_{n_cigar} { }

  enum class RegionType : unsigned char {
    BEGIN,
    ON_GENOME,
    OFF_GENOME_JUNCTION,
    OFF_GENOME_OTHER,
    END
  };

  struct Region {
    ClosedInterval coordinates_;
    uint32_t position_;  // alignment offset
    RegionType type_;

    Region() : coordinates_{}, position_{}, type_{RegionType::END} { }
    Region(ClosedInterval coordinates, uint32_t position, RegionType type)
        : coordinates_{coordinates}, position_{position}, type_{type} { }
    explicit Region(position_t begin_at)
        : coordinates_{ClosedInterval::FromStartLength(begin_at, 0)},
          position_{0},
          type_{RegionType::END} { }
    Region(const Region&) = default;
    Region(Region&&) = default;
    Region& operator=(const Region&) = default;
    Region& operator=(Region&&) = default;
  };

  class Iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = Region;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

   private:
    const CigarRegions& parent_;
    uint32_t idx_;
    value_type region_;
    uint32_t prev_dposition_;  // previous change in alignment offset

    // apply next cigar operation (increment_idx = 0 used for first index)
    template <uint32_t increment_idx>
    void apply_next_op() {
      idx_ = std::min(parent_.n_cigar_, idx_ + increment_idx);
      if (idx_ == parent_.n_cigar_) {
        // we don't use END type -- we just use idx_, no need to do this
        // region_.type_ = RegionType::END;
        // also could update position/coordinates but would again be pointless
        return;
      }
      const char cigar_op = bam_cigar_op(parent_.cigar_[idx_]);
      const char cigar_type = bam_cigar_type(cigar_op);
      if (cigar_type == 0) {
        // does not advance either reference/query so we can just ignore
        apply_next_op<1>();
        return;
      }
      const uint32_t cigar_oplen = bam_cigar_oplen(parent_.cigar_[idx_]);
      // update position, next change in position
      region_.position_ += prev_dposition_;
      prev_dposition_
        = cigar_type & CIGAR_CONSUMES_QUERY ? cigar_oplen : 0;
      // update coordinates
      region_.coordinates_ = ClosedInterval::FromStartLength(
          1 + region_.coordinates_.end,
          cigar_type & CIGAR_CONSUMES_REFERENCE ? cigar_oplen : 0);
      // update type for region
      region_.type_
        = cigar_type == CIGAR_CONSUMES_REFERENCE
          ? (cigar_op == BAM_CREF_SKIP
              ? RegionType::OFF_GENOME_JUNCTION
              : RegionType::OFF_GENOME_OTHER)
          : RegionType::ON_GENOME;
    }

   public:
    Iterator(const CigarRegions& parent, bool begin)
        : parent_{parent},
          idx_{begin ? 0 : parent_.n_cigar_},
          region_{begin ? Region{1 + parent_.genomic_pos_} : Region{}},
          prev_dposition_{0} {
      apply_next_op<0>();
    }
    Iterator& operator++() {
      apply_next_op<1>();
      return *this;
    }

    const reference operator*() noexcept { return region_; }
    pointer operator->() noexcept { return &region_; }

    bool operator==(const Iterator& rhs) const {
      return idx_ == rhs.idx_;
    }
    bool operator!=(const Iterator& rhs) const {
      return !(*this == rhs);
    }
  };

  Iterator begin() { return Iterator{*this, true}; }
  Iterator end() { return Iterator{*this, false}; }
};

}  // namespace bam
}  // namespace majiq


#endif  // MAJIQ_BAM_CIGARREGIONS_HPP
