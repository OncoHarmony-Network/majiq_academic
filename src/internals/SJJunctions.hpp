/**
 * SJJunctions.hpp
 *
 * Contig regions with different data to represent junctions in different
 * settings
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_SJJUNCTIONS_HPP
#define MAJIQ_SJJUNCTIONS_HPP

#include <vector>
#include <map>
#include <utility>
#include <memory>
#include <stdexcept>
#include <algorithm>

#include "MajiqTypes.hpp"
#include "ContigRegion.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"

namespace majiq {

// did a junction pass different thresholds?
enum class JunctionPassedStatus : unsigned char {
  NOT_PASSED,
  ANNOTATED_PASSED,
  DENOVO_PASSED
};

// total reads and nonzero positions for a junction
struct ExperimentCounts {
  junction_ct_t numreads;
  junction_pos_t numpos;
  ExperimentCounts() : numreads{}, numpos{} { }
  ExperimentCounts(junction_ct_t _numreads, junction_pos_t _numpos)
      : numreads{_numreads}, numpos{_numpos} { }
  ExperimentCounts(const ExperimentCounts&) = default;
  ExperimentCounts(ExperimentCounts&&) = default;
  ExperimentCounts& operator=(const ExperimentCounts&) = default;
  ExperimentCounts& operator=(ExperimentCounts&&) = default;
  JunctionPassedStatus passed(junction_ct_t minreads, junction_ct_t mindenovo,
      junction_pos_t minpos) const noexcept {
    // NOTE: assumes mindenovo >= minreads
    if (numpos < minpos || numreads < minreads) {
      return JunctionPassedStatus::NOT_PASSED;
    } else if (numreads >= mindenovo) {
      return JunctionPassedStatus::DENOVO_PASSED;
    } else {  // minreads <= numreads < mindenovo
      return JunctionPassedStatus::ANNOTATED_PASSED;
    }
  }
};


class SJJunction : public detail::ContigRegion<OpenInterval, ExperimentCounts> {
  using BaseT = detail::ContigRegion<OpenInterval, ExperimentCounts>;

 public:
  const junction_ct_t& numreads() const noexcept { return data.numreads; }
  junction_ct_t& numreads() noexcept { return data.numreads; }
  const junction_pos_t& numpos() const noexcept { return data.numpos; }
  junction_pos_t& numpos() noexcept { return data.numpos; }
  JunctionPassedStatus passed(junction_ct_t minreads, junction_ct_t mindenovo,
      junction_pos_t minpos) const noexcept {
    return data.passed(minreads, mindenovo, minpos);
  }

  SJJunction(KnownContig contig, OpenInterval coordinates,
      GeneStrandness strand, ExperimentCounts counts)
      : BaseT{contig, coordinates, strand, counts} { }
  SJJunction(KnownContig contig, OpenInterval coordinates,
      GeneStrandness strand)
      : SJJunction{contig, coordinates, strand, ExperimentCounts{}} { }
  SJJunction() : SJJunction{KnownContig{}, OpenInterval{}, GeneStrandness{}} { }
  SJJunction(const SJJunction&) = default;
  SJJunction(SJJunction&&) = default;
  SJJunction& operator=(const SJJunction&) = default;
  SJJunction& operator=(SJJunction&&) = default;
};
using SJJunctions = detail::Regions<SJJunction, true>;

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONS_HPP
