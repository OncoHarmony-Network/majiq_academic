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

#include <optional>
#include <algorithm>

#include "../MajiqTypes.hpp"
#include "../Interval.hpp"
#include "CigarUtils.hpp"

namespace majiq {
namespace bam {

class CigarRegions : public CigarBase {
 public:
  CigarRegions(
      position_t genomic_pos, int32_t read_length,
      uint32_t* cigar, uint32_t n_cigar)
      : CigarBase(genomic_pos, read_length, cigar, n_cigar) { }

  enum class RegionType : unsigned char {
    BEGIN,
    ON_GENOME,
    OFF_GENOME_JUNCTION,
    OFF_GENOME_OTHER,
    END
  };

  struct Region {
    ClosedInterval coordinates_;
    position_t position_;  // alignment offset
    RegionType type_;

    Region() : coordinates_{}, position_{-1}, type_{RegionType::END} { }
    Region(ClosedInterval coordinates, position_t position, RegionType type)
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
        = cigar_type & detail::CIGAR_CONSUMES_QUERY ? cigar_oplen : 0;
      // update coordinates
      region_.coordinates_ = ClosedInterval::FromStartLength(
          1 + region_.coordinates_.end,
          cigar_type & detail::CIGAR_CONSUMES_REFERENCE ? cigar_oplen : 0);
      // update type for region
      region_.type_
        = cigar_type == detail::CIGAR_CONSUMES_REFERENCE
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
    const reference operator*() noexcept {
      return region_;
    }
    Iterator& operator++() {
      apply_next_op<1>();
      return *this;
    }
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
