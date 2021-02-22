/**
 * SJJunctionsPositions.hpp
 *
 * Underlying raw per-position coverage for SJJunctions
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_SJJUNCTIONSPOSITIONS_HPP
#define MAJIQ_SJJUNCTIONSPOSITIONS_HPP

#include <vector>
#include <memory>
#include <tuple>
#include <stdexcept>

#include "SJJunctions.hpp"
#include "Exons.hpp"
#include "GeneIntrons.hpp"
#include "ContigIntrons.hpp"
#include "MajiqConstants.hpp"


namespace majiq {
// identity/read count for nonzero junction positions
struct PositionReads {
  junction_pos_t pos;
  junction_ct_t reads;
};
inline bool operator<(const PositionReads& x, const PositionReads& y) noexcept {
  // sorted by reads first
  return std::tie(x.reads, x.pos) < std::tie(y.reads, y.pos);
}

class SJJunctionsPositions {
 private:
  std::shared_ptr<SJJunctions> junctions_;
  std::vector<PositionReads> reads_;
  std::vector<size_t> offsets_;
  const junction_pos_t num_positions_;

  void check_valid() const;

 public:
  static SJJunctionsPositions FromBam(
      const char* infile, ExperimentStrandness exp_strandness, int nthreads);

  size_t num_junctions() const noexcept { return junctions_->size(); }
  size_t size() const noexcept { return reads_.size(); }
  junction_pos_t num_positions() const { return num_positions_; }
  const std::shared_ptr<SJJunctions>& junctions() { return junctions_; }
  const std::vector<PositionReads>& reads() { return reads_; }
  const std::vector<size_t>& offsets() { return offsets_; }

  SJJunctionsPositions()
      : junctions_{std::make_shared<SJJunctions>()},
        reads_{},
        offsets_{0},
        num_positions_{0} { }
  SJJunctionsPositions(
      const std::shared_ptr<SJJunctions>& junctions,
      const std::vector<PositionReads>& reads,
      const std::vector<size_t>& offsets,
      junction_pos_t num_positions)
      : junctions_{junctions},
        reads_{reads},
        offsets_{offsets},
        num_positions_{num_positions} { check_valid(); }
  SJJunctionsPositions(
      std::shared_ptr<SJJunctions>&& junctions,
      std::vector<PositionReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t num_positions)
      : junctions_{junctions},
        reads_{reads},
        offsets_{offsets},
        num_positions_{num_positions} { check_valid(); }
  SJJunctionsPositions(const SJJunctionsPositions& x) = default;
  SJJunctionsPositions(SJJunctionsPositions&& x) = default;
  SJJunctionsPositions& operator=(const SJJunctionsPositions& x) = delete;
  SJJunctionsPositions& operator=(SJJunctionsPositions&& x) = delete;
};

class SJIntronsBins {
 private:
  const ContigIntrons introns_;
  const std::vector<PositionReads> reads_;
  const std::vector<size_t> offsets_;
  const junction_pos_t num_bins_;

  void check_valid() const;  // raise error if defined inconsistently

 public:
  static SJIntronsBins FromBam(
      const char* infile, const junction_pos_t num_bins, const Exons& exons,
      const GeneIntrons& gene_introns, ExperimentStrandness exp_strandness,
      int nthreads);

  size_t num_introns() const noexcept { return introns_.size(); }
  size_t size() const noexcept { return reads_.size(); }
  junction_pos_t num_bins() const noexcept { return num_bins_; }
  const ContigIntrons& introns() { return introns_; }
  const std::vector<PositionReads>& reads() { return reads_; }
  const std::vector<size_t>& offsets() { return offsets_; }

  SJIntronsBins(
      const ContigIntrons& introns,
      const std::vector<PositionReads>& reads,
      const std::vector<size_t>& offsets,
      junction_pos_t num_bins)
      : introns_{introns},
        reads_{reads},
        offsets_{offsets},
        num_bins_{num_bins} {
    check_valid();
  }
  SJIntronsBins(
      ContigIntrons&& introns,
      std::vector<PositionReads>&& reads,
      std::vector<size_t>&& offsets,
      junction_pos_t num_bins)
      : introns_{introns},
        reads_{reads},
        offsets_{offsets},
        num_bins_{num_bins} {
    check_valid();
  }
  SJIntronsBins(const SJIntronsBins&) = default;
  SJIntronsBins(SJIntronsBins&&) = default;
  SJIntronsBins& operator=(const SJIntronsBins&) = delete;
  SJIntronsBins& operator=(SJIntronsBins&&) = delete;
};

}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONSPOSITIONS_HPP
