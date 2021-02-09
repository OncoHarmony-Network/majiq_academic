/**
 * CigarJunctions
 *
 * Class for extracting junctions from CIGAR string
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_BAM_CIGARJUNCTIONS_HPP
#define MAJIQ_BAM_CIGARJUNCTIONS_HPP

#include <htslib/sam.h>
#include <utility>
#include "../MajiqTypes.hpp"
#include "../Interval.hpp"
#include "CigarUtils.hpp"

namespace majiq {
namespace bam {

template<uint32_t min_overhang>
class CigarJunctions : public CigarBase {
 public:
  CigarJunctions(
      position_t genomic_pos, int32_t read_length,
      uint32_t* cigar, uint32_t n_cigar)
      : CigarBase(genomic_pos, read_length, cigar, n_cigar) {
  }

  class Iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = std::pair<OpenInterval, junction_pos_t>;
    using difference_type = int;
    using pointer = value_type*;
    using reference = value_type&;

   private:
    // underlying CigarJunctions
    const CigarJunctions<min_overhang>& x_;
    // current offsets
    uint32_t cigar_ndx_;
    uint32_t genomic_offset_;
    uint32_t alignment_offset_;
    // current junction coordinates, if any, and position
    OpenInterval junction_;
    junction_pos_t position_;

    // current query offset (on read) includes left overhang
    uint32_t read_offset() const noexcept {
      return alignment_offset_ + x_.clipping_lengths_.first;
    }

    // get junction associated with current cigar position/offsets
    void UpdateJunctionInformation() {
      junction_.start
        = x_.genomic_pos_ + genomic_offset_;
      junction_.end
        = 1 + junction_.start + bam_cigar_oplen(x_.cigar_[cigar_ndx_]);
      position_ = read_offset() - min_overhang;
    }

   public:
    const value_type operator*() const noexcept {
      return std::make_pair(junction_, position_);
    }
    bool operator==(const Iterator& rhs) noexcept {
      return cigar_ndx_ == rhs.cigar_ndx_ && x_.cigar_ == rhs.x_.cigar_;
    }
    bool operator!=(const Iterator& rhs) noexcept {
      return !(*this == rhs);
    }

    /**
     * Move iterator to next junction or end of cigar string
     * @note should UpdateJunctionInformation when ends at junction
     */
    Iterator& operator++() {
      // assume that we've already visited a junction at cigar_ndx_ only
      // violated when min_overhang = 0 in nonrealistic case that is handled
      // by the constructor
      bool first = true;
      // update offsets and move cigar_ndx_ forward until next junction
      for (; cigar_ndx_ < x_.n_cigar_; ++cigar_ndx_) {
        const char cigar_op = bam_cigar_op(x_.cigar_[cigar_ndx_]);
        const uint32_t cigar_oplen = bam_cigar_oplen(x_.cigar_[cigar_ndx_]);
        // currently (1) a junction, (2) not first iteration, (3) past overhang
        if (cigar_op == BAM_CREF_SKIP
            && !first
            && alignment_offset_ >= min_overhang) {
          UpdateJunctionInformation();  // update junction_
          break;
        } else {
          first = false;
        }
        // cigar_type: bit 1 if consumes query, 2 if consumes reference
        const char cigar_type = bam_cigar_type(cigar_op);
        // update offset on query (i.e. alignment)
        if (cigar_type & detail::CIGAR_CONSUMES_QUERY) {
          alignment_offset_ += cigar_oplen;
          if (alignment_offset_ > x_.alignment_length_ - min_overhang) {
            // future operations can not lead to valid junction
            // set cigar_ndx_ to end and quit
            cigar_ndx_ = x_.n_cigar_;
            break;
          }
        }
        // update offset on reference (i.e. genome)
        if (cigar_type & detail::CIGAR_CONSUMES_REFERENCE) {
          genomic_offset_ += cigar_oplen;
        }
      }
      // when done, return reference to self
      return *this;
    }

    Iterator(const CigarJunctions<min_overhang>& x, bool begin)
        : x_{x},
          cigar_ndx_{begin ? 0 : x_.n_cigar_},
          genomic_offset_{0},
          alignment_offset_{0},
          junction_{},
          position_{0} {
      // if end don't need to do anything, but otherwise get first junction
      if (begin) {
        if constexpr(min_overhang <= 0) {
          // theoretical special case if we don't have overhang: junction on
          // first cigar operation. operator++ assumes junction on current index
          // has already been seen, so we need to skip it
          const char first_op = bam_cigar_op(x_.cigar_[cigar_ndx_]);
          if (first_op == BAM_CREF_SKIP) {
            UpdateJunctionInformation();  // update using current junction
          }
        }
        // otherwise, just use our increment operator
        ++(*this);
      }
    }
  };

  // forward declaration of iterator start/end over cigar junctions
  Iterator begin() { return Iterator(*this, true); }
  Iterator end() { return Iterator(*this, false); }
};

}  // namespace bam
}  // namespace majiq

#endif  // MAJIQ_BAM_CIGARJUNCTIONS_HPP
