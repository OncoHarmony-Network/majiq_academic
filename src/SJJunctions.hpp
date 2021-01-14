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
#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <string>

#include "MajiqTypes.hpp"
#include "Regions.hpp"
#include "Interval.hpp"
#include "Contigs.hpp"
#include "bam/Alignments.hpp"
#include "bam/CigarJunctions.hpp"

namespace majiq {
// total reads and nonzero positions for a junction
struct ExperimentCounts {
  junction_ct_t numreads;
  junction_pos_t numpos;
};
// group passed (build, denovo)
struct GroupPassed {
  bool build_passed;
  bool denovo_passed;
};


using ContigJunction = detail::ContigRegion<OpenInterval, detail::EmptyDataT>;
using SJJunction = detail::ContigRegion<OpenInterval, ExperimentCounts>;
using GroupJunction = detail::ContigRegion<OpenInterval, GroupPassed>;

using SJJunctions = detail::ContigRegions<SJJunction>;
using GroupJunctions = detail::ContigRegions<GroupJunction>;

// how do we get from SJJunctions to GroupJunctions?
class GroupJunctionsGenerator {
 private:
  const junction_ct_t minreads_;
  const junction_ct_t mindenovo_;  // REQUIRED >= minreads
  const junction_pos_t minpos_;

  std::shared_ptr<Contigs> contigs_;

  // number of experiments passing filters
  struct GroupCounts {
    uint32_t numpassed;  // passing build filters
    uint32_t numdenovo;  // passing denovo filters
  };
  size_t num_experiments_;
  std::vector<  // over contigs
    std::map<std::pair<OpenInterval, GeneStrandness>, GroupCounts>> counts_;

 public:
  GroupJunctionsGenerator(junction_ct_t minreads, junction_ct_t mindenovo,
      junction_pos_t minpos, std::shared_ptr<Contigs> contigs)
      : minreads_{minreads}, mindenovo_{mindenovo}, minpos_{minpos},
        contigs_{contigs == nullptr ? std::make_shared<Contigs>() : contigs} {
    if (minreads > mindenovo) {
      throw std::invalid_argument("mindenovo must be at least minreads");
    }
  }
  GroupJunctionsGenerator(junction_ct_t minreads, junction_ct_t mindenovo,
      junction_pos_t minpos)
      : GroupJunctionsGenerator{minreads, mindenovo, minpos, nullptr} {
  }
  GroupJunctionsGenerator(const GroupJunctionsGenerator& x) = default;
  GroupJunctionsGenerator(GroupJunctionsGenerator&& x) = default;
  GroupJunctionsGenerator& operator=(
      const GroupJunctionsGenerator& x) = delete;
  GroupJunctionsGenerator& operator=(GroupJunctionsGenerator&& x) = delete;

  /**
   * Create group junctions with data aggregated from input experiments,
   * filtering unpassed junctions using min_experiment
   */
  GroupJunctions PassedJunctions(size_t min_experiments);
  /**
   * Add junctions from experiment to group towards group filters
   */
  void AddExperiment(const SJJunctions& experiment);
};

// identity/read count for nonzero junction positions
struct PositionReads {
  junction_pos_t pos;
  junction_ct_t reads;
};
inline bool operator<(const PositionReads& x, const PositionReads& y) noexcept {
  // sorted by reads first
  return std::tie(x.reads, x.pos) < std::tie(y.reads, y.pos);
}

class SJPositionReads {
 private:
  std::vector<PositionReads> reads_;
  std::vector<size_t> offsets_;

  static bool is_valid(const std::vector<PositionReads>& reads,
      const std::vector<size_t>& offsets) {
    if (offsets.empty() || offsets.back() != reads.size()) { return false; }
    return offsets.end() == std::adjacent_find(offsets.begin(), offsets.end(),
        [&reads](size_t start, size_t end) {
        // found invalid if offsets decreasing or reads not sorted
        return end < start
          || !std::is_sorted(reads.begin() + start, reads.begin() + end);
        });
  }

 public:
  inline size_t num_junctions() const noexcept { return offsets_.size() - 1; }
  const std::vector<PositionReads>& reads() { return reads_; }
  const std::vector<size_t>& offsets() { return offsets_; }

  SJPositionReads() : reads_{}, offsets_{0} { }
  SJPositionReads(const std::vector<PositionReads>& reads,
      const std::vector<size_t>& offsets) {
    if (!is_valid(reads, offsets)) {
      throw std::invalid_argument("SJPositionReads given invalid arguments");
    }
    reads_ = reads;
    offsets_ = offsets;
  }
  SJPositionReads(std::vector<PositionReads>&& reads,
      std::vector<size_t>&& offsets) {
    if (!is_valid(reads, offsets)) {
      throw std::invalid_argument("SJPositionReads given invalid arguments");
    }
    reads_ = reads;
    offsets_ = offsets;
  }
  SJPositionReads(const SJPositionReads& x) = default;
  SJPositionReads(SJPositionReads&& x) = default;
  SJPositionReads& operator=(const SJPositionReads& x) = default;
  SJPositionReads& operator=(SJPositionReads&& x) = default;
};

class SJJunctionsAll {
 private:
  std::shared_ptr<SJJunctions> junctions_;
  std::shared_ptr<SJPositionReads> reads_;
  junction_pos_t num_positions_;

  static bool is_valid(const std::shared_ptr<SJJunctions>& junctions,
      const std::shared_ptr<SJPositionReads>& reads) {
    return junctions == nullptr
      || reads == nullptr
      || junctions->size() == reads->num_junctions();
  }

 public:
  void set_junctions(std::shared_ptr<SJJunctions> junctions) {
    if (!is_valid(junctions, reads_)) {
      throw std::runtime_error("SJJunctionsAll mismatch of junctions/reads");
    }
    junctions_ = junctions;
  }
  void set_reads(std::shared_ptr<SJPositionReads> reads) {
    if (!is_valid(junctions_, reads)) {
      throw std::runtime_error("SJJunctionsAll mismatch of junctions/reads");
    }
    reads_ = reads;
  }
  junction_pos_t num_positions() const { return num_positions_; }
  const std::shared_ptr<SJJunctions>& junctions() { return junctions_; }
  const std::shared_ptr<SJPositionReads>& reads() { return reads_; }

  SJJunctionsAll(
      const std::shared_ptr<SJJunctions>& junctions,
      const std::shared_ptr<SJPositionReads>& reads,
      junction_pos_t num_positions) {
    if (!is_valid(junctions, reads)) {
      throw std::runtime_error("SJJunctionsAll mismatch of junctions/reads");
    }
    junctions_ = junctions;
    reads_ = reads;
    // NOTE not checking that num_positions is compatible right now
    num_positions_ = num_positions;
  }
  SJJunctionsAll() = default;
  SJJunctionsAll(const SJJunctionsAll& x) = default;
  SJJunctionsAll(SJJunctionsAll&& x) = default;
  SJJunctionsAll& operator=(const SJJunctionsAll& x) = default;
  SJJunctionsAll& operator=(SJJunctionsAll&& x) = default;
};


template <uint32_t min_overhang>
SJJunctionsAll JunctionsFromBam(
    const char* infile, int nthreads, ExperimentStrandness exp_strandness) {
  auto contigs = std::make_shared<Contigs>();
  junction_pos_t max_read_length = 0;
  size_t num_junctions = 0;
  size_t num_junction_positions = 0;
  size_t total_junction_positions = 0;  // how many
  // counts per contig, junction, strand, position
  std::vector<  // per contig
    std::map<  // per junction/strand
      std::pair<OpenInterval, GeneStrandness>,
      std::map<junction_pos_t, junction_ct_t>
    >
  > counts;
  // open up bam and process
  {
    // open up bam for reading with specified number of threads
    bam::AlignmentFile in(infile, nthreads);

    // track contigs -- just copy over contigs from alignment file
    std::for_each(in.target_name(), in.target_name() + in.n_targets(),
        [&contigs](char* seqid) { contigs->add(std::string{seqid}); });
    // expand junctions for all contigs
    counts.resize(in.n_targets());

    // iterate over alignments
    bam::AlignmentRecord aln{};
    int r = 0;
    // while we successfully read in a new alignment
    while ((r = in.NextRead(aln)) >= 0) {
      // update max read length
      max_read_length = std::max(
          max_read_length, static_cast<junction_pos_t>(aln.read_length()));
      // only process if unique/mapped
      if (!aln.unique_mapped()) { continue; }
      // parse junctions from uniquely mapped alignments
      // process cigar for junctions
      auto junctions = aln.cigar_junctions<min_overhang>();
      auto cur_junc = junctions.begin();
      // only process if there are any valid junctions (given overhang)
      if (cur_junc == junctions.end()) { continue; }
      // at least one junction
      auto& contig_counts = counts[aln.tid()];  // index into contig
      GeneStrandness strand = aln.strandness(exp_strandness);
      // add each junction in this alignment
      for (; cur_junc != junctions.end(); ++cur_junc) {
        auto& [interval, position] = *cur_junc;
        // counts for junction
        auto& junction_counts = contig_counts[
          std::make_pair(interval, strand)];
        // increment num_junctions if we haven't seen it before
        if (junction_counts.empty()) { ++num_junctions; }
        // increment position, and num_junction_positions if new
        if (1 == ++junction_counts[position]) { ++num_junction_positions; }
      }  // added each junction
    }  // while successfully loading next read
    if (r < -1) {
      std::ostringstream oss;
      oss << "Error parsing read from BAM " << infile
        << " (error code: " << r << ")";
      // throw std::runtime_error(oss.str());
      std::cerr << oss.str() << std::endl;
    }  // otherwise r = 0 and at end of file
  }  // done processing BAM file
  // number of valid positions (note that type may be unsigned)
  junction_pos_t eff_len
    = (max_read_length + 1 > 2 * min_overhang)
    ? max_read_length + 1 - 2 * min_overhang : 0;
  // initialize vectors for outputs
  std::vector<SJJunction> junctions(num_junctions);
  std::vector<PositionReads> reads(num_junction_positions);
  std::vector<size_t> reads_offsets(num_junctions + 1);
  size_t jidx = 0;
  size_t jpidx = 0;
  for (size_t contig_idx = 0; contig_idx < counts.size(); ++contig_idx) {
    if (counts[contig_idx].empty()) { continue; }
    KnownContig contig = (*contigs)[contig_idx];
    for (auto& [pos_strand, position_reads] : counts[contig_idx]) {
      const junction_pos_t num_positions = position_reads.size();
      // copy over position reads and sort in increasing order of read counts
      std::transform(position_reads.begin(), position_reads.end(),
          reads.begin() + jpidx,
          [](const std::pair<junction_pos_t, junction_ct_t>& x) {
          return PositionReads{x.first, x.second};
          });
      std::sort(reads.begin() + jpidx, reads.begin() + jpidx + num_positions);
      // summarize total reads
      const junction_ct_t num_reads = std::accumulate(
          reads.begin() + jpidx, reads.begin() + jpidx + num_positions,
          junction_ct_t{},
          [](junction_ct_t acc, const PositionReads& x) {
          return acc + x.reads;
          });
      // update junctions, offsets
      junctions[jidx] = SJJunction{contig, pos_strand.first, pos_strand.second,
                                   ExperimentCounts{num_reads, num_positions}};
      reads_offsets[jidx + 1] = reads_offsets[jidx] + num_positions;
      // increment jidx, jpidx
      ++jidx;
      jpidx += num_positions;
    }
  }
  // finally, create desired result
  return SJJunctionsAll{
    std::make_shared<SJJunctions>(junctions),
    std::make_shared<SJPositionReads>(reads, reads_offsets),
    eff_len
  };
}

GroupJunctions GroupJunctionsGenerator::PassedJunctions(
    size_t min_experiments) {
  std::vector<GroupJunction> passed;
  for (size_t contig_idx = 0; contig_idx < contigs_->size(); ++contig_idx) {
    for (auto&& [coord_strand, ct] : counts_[contig_idx]) {
      if (ct.numpassed >= min_experiments) {
        passed.emplace_back(
            (*contigs_)[contig_idx], coord_strand.first, coord_strand.second,
            GroupPassed{ct.numpassed >= min_experiments,
                        ct.numdenovo >= min_experiments});
      }
    }
  }
  // create GroupJunctions
  return GroupJunctions{passed};
}
void GroupJunctionsGenerator::AddExperiment(const SJJunctions& experiment) {
  for (const SJJunction& junction : experiment) {
    // need min positions to pass no matter what
    if (junction.data.numpos < minpos_) { continue; }
    // otherwise passed in general vs for denovo
    const bool minreads_passed = junction.data.numreads >= minreads_;
    const bool mindenovo_passed = junction.data.numreads >= mindenovo_;
    if (mindenovo_passed || minreads_passed) {
      // need to update counts_
      size_t contig_idx = contigs_->add(junction.contig.get());
      if (contig_idx >= counts_.size()) {
        // make sure index available for contig
        counts_.resize(contig_idx + 1);
      }
      // get reference to value for contig/coordinates/strand to add
      // make two copies of unstranded junctions to handle experiments with
      // mixed strandedness (unstranded/stranded)
      std::vector<GeneStrandness> strands;
      if (junction.strand == GeneStrandness::AMBIGUOUS) {
        // both strands in sorted order
        if constexpr(GeneStrandness::FORWARD < GeneStrandness::REVERSE) {
          strands = {GeneStrandness::FORWARD, GeneStrandness::REVERSE};
        } else {
          strands = {GeneStrandness::REVERSE, GeneStrandness::FORWARD};
        }
      } else {
        strands = {junction.strand};
      }
      for (GeneStrandness strand : strands) {
        auto& counts = counts_[contig_idx][
          std::make_pair(junction.coordinates, strand)];
        if (minreads_passed) {
          ++(counts.numpassed);
        }
        if (mindenovo_passed) {
          ++(counts.numdenovo);
        }
      }
    }  // done adding junction that passed
  }  // end loop over junctions
  // update num_experiments_
  ++num_experiments_;
}
}  // namespace majiq

#endif  // MAJIQ_SJJUNCTIONS_HPP
