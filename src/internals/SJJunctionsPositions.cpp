/**
 * SJJunctionsPositions.cpp
 *
 * Implementation of non-inline functions in SJJunctionsPositions
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "SJJunctionsPositions.hpp"

#include <memory>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <utility>

#include "SJJunctions.hpp"
#include "bam/Alignments.hpp"
#include "bam/CigarJunctions.hpp"


namespace majiq {

bool SJJunctionsPositions::is_valid(
    const std::shared_ptr<SJJunctions>& junctions,
    const std::vector<PositionReads>& reads,
    const std::vector<size_t>& offsets,
    junction_pos_t num_positions) {
  if (
      // must have junctions
      junctions == nullptr
      // offsets[i], offsets[i + 1] correspond to junctions[i] --> match sizes
      || junctions->size() != offsets.size() - 1
      // last offset must correspond to size of reads
      || offsets.back() != reads.size()) {
    return false;
  }
  // for each junction
  for (size_t jidx = 0; jidx < junctions->size(); ++jidx) {
    // offsets into junction positions
    auto jp_start = reads.begin() + offsets[jidx];
    auto jp_end = reads.begin() + offsets[jidx + 1];
    if (
        // offsets must be non-decreasing
        jp_end < jp_start
        // difference in offsets correspond to numpos
        || jp_end - jp_start != (*junctions)[jidx].data.numpos
        // can't have numpos greater than total possible
        || (*junctions)[jidx].data.numpos > num_positions
        // all junction positions must be in [0, num_positions)
        || !std::all_of(jp_start, jp_end,
          [num_positions](const PositionReads& x) {
          return x.pos < num_positions;
          })
        // sum of position reads equals junction numreads
        || std::accumulate(jp_start, jp_end, junction_ct_t{},
          [](junction_ct_t s, const PositionReads& x) { return s + x.reads; })
          != (*junctions)[jidx].data.numreads
        // PositionReads in sorted order (by number of reads)
        || !std::is_sorted(jp_start, jp_end)) {
      return false;
    }
  }
  return true;
}

SJJunctionsPositions SJJunctionsPositions::FromBam(
    const char* infile, ExperimentStrandness exp_strandness, int nthreads) {
  auto contigs = std::make_shared<Contigs>();
  junction_pos_t max_read_length = 0;
  size_t num_junctions = 0;
  size_t num_junction_positions = 0;
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
        [&contigs](char* seqid) { contigs->add(Contig{std::string{seqid}}); });
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
      auto junctions = aln.cigar_junctions<USE_MIN_OVERHANG>();
      auto cur_junc = junctions.begin();
      // only process if there are any valid junctions (given overhang)
      if (cur_junc == junctions.end()) { continue; }
      // at least one junction
      auto& contig_counts = counts[aln.tid()];  // index into contig
      GeneStrandness strand = aln.strandness(exp_strandness);
      // add each junction in this alignment
      for (; cur_junc != junctions.end(); ++cur_junc) {
        auto& [coordinates, position] = *cur_junc;
        // counts for junction
        auto& junction_counts = contig_counts[
          std::make_pair(coordinates, strand)];
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
    = (max_read_length + 1 > 2 * USE_MIN_OVERHANG)
    ? max_read_length + 1 - 2 * USE_MIN_OVERHANG : 0;
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
  return SJJunctionsPositions{
    std::make_shared<SJJunctions>(std::move(junctions)), reads, reads_offsets, eff_len
  };
}

}  // namespace majiq
