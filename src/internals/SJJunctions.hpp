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
#include "ExperimentThresholds.hpp"

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
  JunctionPassedStatus passed(
      const ExperimentThresholds& thresholds) const noexcept {
    if (numpos < thresholds.minpos_ || numreads < thresholds.minreads_) {
      return JunctionPassedStatus::NOT_PASSED;
    } else if (numreads >= thresholds.mindenovo_) {
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
  JunctionPassedStatus passed(
      const ExperimentThresholds& thresholds) const noexcept {
    return data.passed(thresholds);
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

class SJJunctions : public detail::Regions<SJJunction, true> {
  using BaseT = detail::Regions<SJJunction, true>;

 public:
  SJJunctions(
      const std::shared_ptr<Contigs>& contigs, std::vector<SJJunction>&& x)
      : BaseT{contigs, std::move(x)} {
    if (parents() == nullptr) {
      throw std::invalid_argument("SJJunctions cannot have null genes");
    }
  }
};

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONS_HPP
